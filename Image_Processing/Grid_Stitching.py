#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
import re
import tempfile
import shutil
import glob

from email.MIMEMultipart import MIMEMultipart
from email.MIMEText import MIMEText
from email.Utils import formatdate
import smtplib

from numpy import zeros

import omero.scripts as scripts
import omero.util.script_utils as script_util
from omero.gateway import BlitzGateway
import omero
from omero.rtypes import *
import omero.util.script_utils as script_utils
import omero.util.pixelstypetopython as pixels_type
from omero.util.tiles import *
from omero.model import *
import omero.cli
from omero.rtypes import wrap
from omero.model import DatasetI, ProjectI

from ome_metadata import OMEExporter


ADMIN_EMAIL = 'admin@omerocloud.qbi.uq.edu.au'
input_dir = ''
output_dir = ''


def delete_tmp(tmp_dir):
    """
    Delete the temporary directory
    
    :param tmp_dir: the path of the directory to be deleted
    """
    try:
        for name in glob.glob("%s/*" % tmp_dir):
            os.remove(name)
        os.rmdir(tmp_dir)
    except:
        pass    


def create_containers(conn, parent, child_image):
    """
    Link the new image to a parent.
    """
    
    updateService = conn.getUpdateService()
    parentDataset = parent.getParent()
    parentProject = parentDataset.getParent()
         
    if parentDataset is None:
        print(
            "No dataset created or found for new images."
            " Images will be orphans."
        )
    else:
        dsLink = omero.model.DatasetImageLinkI()
        dsLink.parent = omero.model.DatasetI(
            parentDataset.getId(), False)
        dsLink.child = omero.model.ImageI(child_image.getId(), False)
        updateService.saveObject(dsLink)
#     if parentProject and parentProject.canLink():
#         # and put it in the   current project
#         projectLink = omero.model.ProjectDatasetLinkI()
#         projectLink.parent = omero.model.ProjectI(
#             parentProject.getId(), False)
#         projectLink.child = omero.model.DatasetI(
#             parentDataset.id.val, False)
#         updateService.saveAndReturnObject(projectLink) 


def list_image_names(conn, ids, file_anns):
    """
    Builds a list of the image names
    
    @param conn: The BlitzGateway connection
    @param ids: Python list of image ids
    """
    image_names = []
    for i,image_id in enumerate(ids):
        img = conn.getObject('Image', image_id)
        if not img:
            continue

        ds = img.getParent()
        if ds:
            pr = ds.getParent()
        else:
            pr = None

        image_names.append(
            "[{0}][{1}] Image {2} : {3} : {4}".format(
            pr and pr.getName() or '-',
            ds and ds.getName() or '-',
            image_id, os.path.basename(img.getName()),
            file_anns[i].getFile().getName())
        )

    return image_names


def email_results(conn,params,image_ids,file_anns):
    """
    E-mail the result to the user.

    @param conn: The BlitzGateway connection
    @param params: The script parameters
    @param image_ids: A python list of the new image omero ids
    """
    if not params['Email_Results']:
        return

    image_names = list_image_names(conn, image_ids, file_anns)

    msg = MIMEMultipart()
    msg['From'] = ADMIN_EMAIL
    msg['To'] = params['Email_address']
    msg['Date'] = formatdate(localtime=True)
    msg['Subject'] = '[OMERO Job] Stitching'
    msg.attach(MIMEText("""
New stitching results files created:

Format:
[parent project/datset][daostorm results] image id : image name : result filename

------------------------------------------------------------------------
%s""" % ("\n".join(image_names))))

    smtpObj = smtplib.SMTP('localhost')
    smtpObj.sendmail(ADMIN_EMAIL, [params['Email_address']], msg.as_string())
    smtpObj.quit()


def get_new_image(conn):    
    log = glob.glob(output_dir + '/stdout.txt')
    with open(log[0],'r') as f:
        ids = f.readlines()
        
    image_id = int(ids[0])
    newImg = conn.getObject('Image',image_id)
    return newImg


def do_import(conn, session, filename, dataset=None, project=None):
    """
    Use the omero command line interface to import the image.
    Use the current session.

    :param conn:            the BlitzGateWay connection
    :param session:         dictionary containing the session ID and
                            hostname
    :param filename:        name of the file being imported to omero
    :param dataset:         the omero dataset where we are importing to
    :param project:         the omero project where we are importing to
    """
    user = conn.getUser()
    group = conn.getGroupFromContext()
    
    sessionId = session['ID']
    hostname = session['host']
    cli = omero.cli.CLI()
    cli.loadplugins()
    cli.invoke(
        ["sessions",
         "login",
         "-s",
         "localhost",
         "-k",
         "{}".format(sessionId)],
        strict=True
    )

    import_args = ["import"]
    if dataset is not None:
        dsId = create_containers(conn, dataset, project)
        import_args.extend(["-d", str(dsId)])

    import_args.append(filename)
    import_args.extend([
        "-s","localhost","-u","{}".format(user.getName())
    ])

    # redirect both stderr and stdout to file
    errlog = output_dir + "/stderr.txt"
    import_args.extend(["---errs",errlog])
    outlog = output_dir + "/stdout.txt"
    import_args.extend(["---file",outlog])
    cli.invoke(import_args, strict=True)
    
    # use stdout to get the id of the new image
    newImg = get_new_image(conn)
    return newImg


def run_imagej_script(stitching_args):
    """
    Here we set-up the ImageJ script and run it from the command line
    using the Docker Image fiji/fiji:latest.
    The script text is written to the temp folder that we're running the script in.
    Note that we also create volumes in the Docker image for input images (tiles)
    and output images (fused).
    
    :param image_name:      filename of image being processed
    :param stitching_args:  a list of arguments for stitching provided by the script gui
    """
    # the script is saved as text file
    with open('imagej_stitching_script.txt', 'r') as file:
        stitching_script_str = file.read()

    # so that the stitching args can be inserted
    stitching_script = stitching_script_str.format(*stitching_args)

    # then it is saved remotely
    script = "stitching.py"
    script_path = os.path.join(input_dir, script)

    # write the macro to a known location that we can pass to ImageJ
    f = open(script_path, 'w')
    f.write(stitching_script)
    f.close()

    # Call dockerized Fiji
    cmd = (
        "docker run --rm -v {0}:/fiji/input -v {1}:/fiji/output"
        "fiji/fiji:latest fiji-linux64 --memory=8000m --headless"
        "'/fiji/input/{2}'".format(input_dir,  output_dir,script)
    )
    os.system(cmd)   


def run_stitching(conn, session, stitching_args):
    """
    Launches the Grid Stitching plugin and uplaods results
    
    :param conn:            the BlitzGateWay connection
    :param session:         dictionary containing the session ID and
                            hostname
    :param stitching_args:  list of arguments for stitching provided by
                            the script gui
    """ 
    
    conn.keepAlive()   
    images = glob.glob(input_dir + '/*.tif')
    
    run_imagej_script(stitching_args)
    stitched = glob.glob('%s/*.tif' % output_dir)
    newImg = None
    if stitched:
        newImg = do_import(conn,session,stitched[0])
        
    return newImg


def download_tiles(conn, image, theC, theZ):
    """ 
    Export every plane in the original image, or those
    selected by the user, as OME-TIFF

    :param conn:    BlitzGateway connection
    :param image:   the omero image being downloaded
    :param theC:    the image channel being downloaded
    :parem theZ:    the z-plane being downloaded
    """
    num_tiles = image.getSizeT()
    image_names = []
    for t in range(num_tiles):
        im_name = "tile_T{}.ome.tif".format(t)
        exporter = OMEExporter(
            conn,
            image,
            input_dir,
            im_name,
            theZ=theZ,
            theC=theC,
            theT=t
        )
        exporter.generate()
        image_names.append(im_name)

    return image_names


def run_processing(conn, session, script_params):
    """
    Collects params and starts the processing
    
    @param conn:          the BlitzGateWay connection
    @param script_params: the parameters collected from the script input
    """ 
        
    global input_dir
    global output_dir

    input_dir = tempfile.mkdtemp(prefix='stitching_input')
    output_dir = tempfile.mkdtemp(prefix='stitching_output')
    
    def empty_dir(dir_path):
        for old_file in os.listdir(dir_path):
            file_path = os.path.join(dir_path, old_file)
            os.unlink(file_path)
        
    message = ""
      
    dataType = script_params["Data_Type"]
    
    # Get the images
    objects, logMessage = script_utils.getObjects(conn, script_params)
    message += logMessage
    if not objects:
        return None, message

    # Concatenate images from datasets
    if dataType == 'Image':
        images = objects
    else:
        images = []
        for ds in objects:
            images += ds.listChildren()
    
    image_ids = [i.getId() for i in images]
    if len(image_ids) > 10:
        message = 'Max number of datasets for batch exceeded (10)'
        return message
    
    new_images = []
    new_ids = []
    for image in conn.getObjects("Image",image_ids):
        # remove input and processed images
        empty_dir(input_dir)
        empty_dir(output_dir)
        theC = None    
        if script_params['Single_Channel']:
            theC = script_params['Channel']

        theZ = None
        sizeZ = image.getSizeZ()
        if script_params['Single_Z']:
            theZ = script_params['Z_slice']
            sizeZ = 1
        if script_params['Range_Z']:
            zstart = script_params['Z_start']
            zstop = script_params['Z_stop']
            theZ = range(zstart,zstop+1)
            sizeZ = len(theZ)
        
        # download the image
        image_names = download_tiles(conn,image,theC,theZ)
                
        results_file = str(image.getId()) + "_stitching.txt"        
        
        stitching_args = (
            results_file,
            script_params['grid_x'],
            script_params['grid_y'],
            script_params['tile_overlap'],
            script_params['fusion_method'],
            script_params['regression_threshold'],
            script_params['ave_displacement_threshold'],
            script_params['abs_displacement_threshold'],
            sizeZ
        )
        
        # do the stitching using dockerized fiji
        new_image = run_stitching(conn,  session,stitching_args)
        
        # if the new image is created
        # link it to the original omero image
        if new_image:
            create_containers(conn, image, new_image)
            new_images.append(new_image)
            new_ids.append(new_image.getId())
     
        if len(new_ids) == 0:
            print("No new images created.")
            return       
  
    if new_images:
        if len(new_images) > 1:
            message += "Created %s new images" % len(new_images)
        else:
            message += "Created a new image"
    else:
        message += "No image created"
        
    
    print(script_params['Email_Results'])
    # if new images were created email the omero user
    if script_params['Email_Results'] and new_images:
        email_results(conn, script_params, new_ids)
    
    # remove downloaded image and outputs
    shutil.rmtree(input_dir)
    shutil.rmtree(output_dir)

    robj = (len(new_images) > 0) and new_images[0]._obj or None
    return robj, message

def validate_email(conn, params):
    """
    Checks that a valid email address is present for the user_id

    @param conn: The BlitzGateway connection
    @param params: The script parameters
    """
    userEmail = ''
    if params['Email_address']:
        userEmail = params['Email_address']
    else:
        user = conn.getUser()
        user.getName()
        # Initialises the proxy object for simpleMarshal
        dic = user.simpleMarshal()
        if 'email' in dic and dic['email']:
            userEmail = dic['email']

    params['Email_address'] = userEmail

    # Validate with a regular expression. Not perfect but it will do
    return re.match(
        "^[a-zA-Z0-9._%-]+@[a-zA-Z0-9._%-]+.[a-zA-Z]{2,6}$",
        userEmail
    )

def run_as_script():
    """
    The main entry point of the script, as called by the client via the
    scripting service, passing the required parameters.
    """
                
    dataTypes = [rstring("Dataset"),rstring("Image")]
    fusion_method = [
        rstring('Linear Blending'), rstring('Average'), rstring('Median'),
        rstring('Max. Intensity'), rstring('Min. Intensity'),
        rstring('Intensity of random input tile'),
        rstring('Do not fuse tiles (only write TileConfiguration')
    ]

    client = scripts.client(
        'Grid_Stitching.py',"""Run the "Stitching" FIJI plugin.
MAXIMUM NUMBER OF DATASETS FOR BATCH IS FIVE!""",

    scripts.String(
        "Data_Type", optional=False, grouping="01",
        description="Choose source of images (only Image supported)",
        values=dataTypes, default="Image"
    ),
        
    scripts.List(
        "IDs", optional=False, grouping="02",
        description="IDs of images to be stitched"
    ).ofType(rlong(0)),
                            
    scripts.Bool(
        "Single_Channel", grouping="03",default=False,
        description=(
            "Stitch all channels or a single channel?"
            "Uncheck for all channels"
        )
    ),
                                                        
    scripts.Int("Channel", grouping="03.1",
        description="channel to be stitched"),

    scripts.Bool(
        "Single_Z", grouping="04",default=False,
        description=(
            "Stitch all z slices or a single"
            "slice? Uncheck for all z slices"
        )
    ),
                                                        
    scripts.Int(
        "Z_slice", grouping="04.1", description="z-slice to be stitched"
    ),
                            
    scripts.Bool(
        "Range_Z", grouping="05",default=False,
        description="Stitch a range of slices? Uncheck for all z slices"
    ),
                                                        
    scripts.Int(
        "Z_start", grouping="05.1", description="start at this z-slice"
    ),
                        
    scripts.Int(
        "Z_stop", grouping="05.2", description="stop at this z slice"
    ),
                                                        
    scripts.Int(
        "grid_x", optional=False, grouping="06", default=2,
        description="how many tiles in the x-direction"
    ),

    scripts.Int(
        "grid_y", optional=False, grouping="07", default=2,
        description="how many tiles in the y-direction"
    ),      
                            
    scripts.Int(
        "tile_overlap", optional=False, grouping="08",
        default=20,
        description="percentage overlap between tiles"
    ),  
                            
    scripts.String(
        "fusion_method", optional=False, grouping="09",
        default='Linear Blending',
        description="method used to fuse the tiles",
        values=fusion_method
    ),
                        
    scripts.Float(
        "regression_threshold", optional=False,
        grouping="10",default=0.3,
        description="global optimisation parameter"
    ),  

    scripts.Float(
        "ave_displacement_threshold", optional=False,
        grouping="11",default=2.5,
        description="global optimisation parameter"
    ),

    scripts.Float(
        "abs_displacement_threshold", optional=False, 
        grouping="12",default=3.5,
        description="global optimisation parameter"
    ),                                             
                            
    scripts.Bool(
        "Email_Results", grouping="13", default=False,
        description="E-mail the results"
    ),
                            
    scripts.String(
        "Email_address", grouping="13.1",
        description="Specify e-mail address"
    ),
        
    authors = ["Daniel Matthews", "QBI"],
    institutions = ["University of Queensland"],
    contact = "d.matthews1@uq.edu.au",
    )

    try:
        client.enableKeepAlive(3600)
        # process the list of args above.
        scriptParams = {}
        for key in client.getInputKeys():
            if client.getInput(key):
                scriptParams[key] = client.getInput(key, unwrap=True)
                
        session = {}
        session['ID'] = client.getSessionId()
        session['host'] = client.getProperty('omero.host')
                
        print(scriptParams)
        
        # wrap client to use the Blitz Gateway
        conn = BlitzGateway(client_obj=client)
        if (scriptParams['Email_Results'] and 
            not validate_email(conn, scriptParams)):

            client.setOutput(
                "Message", rstring("No valid email address")
            )
            return
        
        # process images in Datasets
        robj,message = run_processing(conn, session, scriptParams)
        client.setOutput("Message", rstring(message))
        if robj is not None:
            client.setOutput("Result", robject(robj))
    finally:
        client.closeSession()

if __name__ == "__main__":
    run_as_script()
