# Several tools for the MDM reduction pipeline

import astropy.io.ascii as at
import astropy.io.fits as fits
import numpy as np
import subprocess
import pdb

def makelist(fname):
    '''
    This functions makes an ascii text file with a list of all fits files in the current folder. Returns None.

    fname - String; Name of output file.
    '''

    subprocess.call("ls *.fit* > " + fname, shell=True)


def prep(ref_lamp, output_file="to_reduce.lis",
         data_section="[10:4086,1530:2567]", bias_section="[301:364,1:2048]"):
    '''
    This function generates an ascii file with information on each image fits file to process.

    ref_lamp - String; filename of reference lamp image.
    output_file - String; name of out put ascii text file.
    data_section - String; specifies the pixel array that covers the data section.
    bias_section - String; scpecifices the pixel array that covers the bias section (ignore it for OSMOS images).
    '''
    subprocess.call("ls *.fit* > prep.lis", shell=True)

    infile = at.read("prep.lis", data_start=0)
    names = infile.columns[0].data

    with open(output_file,"w") as f:
        for i, name in enumerate(names):
            hdu = fits.open(name)

            # For OSMOS files, skip fits files that are just sky images
            try:
                tmpslit = hdu[0].header["SLITID"]
            except KeyError:
                tmpslit = None
            if tmpslit is not None:
                if tmpslit.strip().upper() == "OPEN": continue

            imgtyp = hdu[0].header["IMAGETYP"].strip().upper()
            #datareg = hdu[0].header["DATASEC"].strip()
            #biasreg = hdu[0].header["BIASSEC"].strip()
            obj_name = hdu[0].header["OBJECT"].strip()
            obj_name = obj_name.replace(" ", "") # No white spaces allowed in obj names

            if imgtyp in ["FOCUS", "LAMP", "LAMPS", "COMP"]:
                category = "lamp"
            elif imgtyp == "BIAS":
                category = "bias"
            elif imgtyp == "OBJECT":
                category = "obj"
            elif "FLAT" in imgtyp:
                category = "flat"
            elif imgtyp == "STANDARD":
                category = "std"
            else:
                print("TYPE {0} NOT KNOWN".format(imgtyp))
                category = ""
            f.write("{:>14} {:>6} {:>18} {:>16} {:>16} {:>18}\n".format(
                    name.rsplit(".", 1)[0], category, obj_name,
                    data_section, bias_section, ref_lamp))
        f.close()

    subprocess.call("rm prep.lis", shell=True)

def read_list(imagelist,
              column_headers=["ccdno","type","target","image_region",
                              "bias_region","reference_lamp"],
              obj_types=["obj","std","flat","bias","lamp"],
              science_types=["obj","std"], return_other_cols=True,
              return_regions=True):
    """
    read in an object list made by prep().

    input
    -----
    imagelist: string, filename
        filename for a list of spectra output by prep()

    column_headers: array_like containing strings

    obj_types: array_like containing strings
        (default = ["obj","std","flat","bias","lamp"])
        which object types to save from everything in imagelist

    science_types: array_like containing strings (default=["obj","std"])
        all the object types that are astronomical objects -
        typically true science targets and flux standards

    return_regions: bool (default=True)
        if True, output will contain strings for the data and overscan regions
        keyed by "good_region" and "overscan_region"

    return_lamps: bool (default=True)
        if True, output will contain a list of reference lamps for the science targets
        keyed by "science_lamps"

    output
    ------
    returns a dictionary keyed by obj_types with the filename bases
    of each requested type, along with
        "science_list" - filenames for all types in science_types
        "science_names" - target names for all types in science_types
        "science_lamps" - reference lamp for all types in science_types
        "std_names" - names associated with "std" (only if "std" is requested)
    if requested, also includes
        "good_region"
        "overscan_region"
    """
    output = {}

    # set up arrays
    for object_type in obj_types:
        output[object_type] = []
        if object_type=="std":
            output["std_names"] = []
    output["science_list"] = []
    output["science_names"] = []
    if return_other_cols and (len(column_headers)>3):
        for col_header in column_headers[3:]:
            output["science_{}".format(col_header)] = []

    # read in list of files
    image_list = at.read(imagelist, converters={'col1':[at.convert_numpy(np.str)]},
                         data_start=0, names=column_headers)

    # save the base filenames for all the requested types, and sometimes names
    for j,target_type in enumerate(image_list["type"]):
        for i,object_type in enumerate(obj_types):
            if target_type == object_type:
                output[object_type].append(image_list["ccdno"][j])
                # if it's a standard, save the standard name for later
                if object_type=="std":
                    output["std_names"].append(image_list["target"][j])
                else:
                    continue
            else:
                continue
        # for final output, standards are "science"
        # so save "obj" and "std" into one list (default types, anyway)
        for k,science_type in enumerate(science_types):
            if target_type == science_type:
                output["science_list"].append(image_list["ccdno"][j])
                output["science_names"].append(image_list["target"][j])
                if return_other_cols:
                    for col_header in column_headers[3:]:
                        output["science_{}".format(col_header)
                               ].append(image_list[col_header][j])
                else:
                    continue
            else:
                continue

    if return_regions:
        output["good_region"] = image_list["image_region"][0]
        output["overscan_region"] = image_list["bias_region"][0]

    return output

def read_reduction_list(imagelist, obj_types=["obj","std","flat","bias","lamp"],
                        science_types=["obj","std"], return_regions=True):
    '''
    Read the list of files to be reduced that was output by prep().
    '''

    output = read_list(imagelist,obj_types=obj_types,science_types=science_types,
                       return_regions=return_regions,return_other_cols=True)
    junk = output.pop("science_image_region")
    junk = output.pop("science_bias_region")

    return output

def read_OI_shifts(shiftlist, science_types=["obj","std"]):
    '''
    Read the list of files that was output by OIshift_corr.py.
    '''
    output = read_list(shiftlist,column_headers=["ccdno","type","target",
                                                 "shift","shift_err","shift_qual"],
                       obj_types=science_types,science_types=science_types,
                       return_regions=False,return_other_cols=True)

    return output

def generate_shift_list(imagelist="to_reduce.lis", output_filename="to_shift.lis", science_types=["obj","std"]):
    '''
    Takes the original list created by prep() and removes all calibration
    files. The output list can be used to run OIshift_corr.pro and fluxcal.py.
    '''

    image_list = at.read(imagelist,data_start=0,
                         names=["ccdno","type","target","image_region",
                                "bias_region","reference_lamp"])

    output_file = open(output_filename,"w")

    for j,target_type in enumerate(image_list["type"]):
        for i,object_type in enumerate(science_types):
            if target_type==object_type:
                output_file.write(" {}  {}  {}\n".format(image_list["ccdno"][j],
                                                       object_type,
                                                       image_list["target"][j]))

    output_file.close()

def check_duplicate_names(science_names):
    '''
    Makes sure that all names in the output target list are unique. Returns a Boolean indicating whether the input list contains duplicated names.

    science_names - List of Strings.
    '''

    duplicates = True

    unique_names,name_counts = np.unique(science_names,return_counts=True)

    if len(unique_names)==len(science_names):
        print("All target names are unique!")
        duplicates = False
    else:
        print("The following names are repeated:")
        print(unique_names[name_counts>1])

    return duplicates

def crop_fits(fname, xsize, ysize, croploc='center', prefix='c_', suffix=None):
    '''
    Crop a fits image using the parameters provided. If file has more than one image, it only considers the first one.

    fname - String, the full path of the fits file; if only a filename is provided, it will look for the file in the current directory.
    xsize - Int, the desired X size (columns) in pixels.
    ysize - Int, the desired Y size (rows) in pixels.
    croploc - ['center'(default), 'centerlow' 'upper right', 'upper left', 'lower left', 'lower right', 'upper center'], set location around which to crop image. If 'center', then it crops image centered in the image center. If 'upper right', then it crops image to size [xsize,ysize] anchored in the upper right corner, and so on. If 'centerlow', it crops the image using the X location of the center as the center of the new image, and the Y location of the center as the lower anchor of the new image.
    prefix - String, prefix to add to new fits file. If both prefix and suffix are None, the original fits file is overwritten with the new one.
    suffix - String, suffix to add to new fits file. If both prefix and suffix are None, the original fits file is overwritten with the new one.

    Returns:
    - the new fits HDU, including the original header information.
    - It also saves a copy of the newly created fits file in the same folder as the original file, with an added suffix to its name, if "suffix" is specified.
    '''

    import os

    # Determine if input is a file or a list of files
    files = []
    if '.fits' in fname:
        files.append(fname)
    else:
        with open(fname, 'r') as infile:
            for line in infile:
                files.append(line[:-1]) # [:-1] is to remove line endings

    for fitsfile in files:
        # Get file path, if provided, and filename
        filepath = fitsfile.rsplit('/',1)[0]
        if filepath == fitsfile:
            filepath = ''
            filename = fitsfile.rsplit('.',1)[0]
        else:
            filepath = filepath + '/'
            filename = fitsfile.rsplit('/',1)[1].rsplit('.',1)[0]

        # Read fits file data
        FitsHDU = fits.open(fitsfile)
        Im = FitsHDU[0].data
        FitsHeader = FitsHDU[0].header
        xsizeorig = FitsHeader['NAXIS1']
        ysizeorig = FitsHeader['NAXIS2']

        # Determine pixel limits for cropping
        if croploc == 'center':
            center = [int(xsizeorig/2), int(ysizeorig/2)]
            xstart = center[0] - int(xsize/2) + 1
            xstop = center[0] + int(xsize/2) + 1
            ystart = center[1] - int(ysize/2)
            ystop = center[1] + int(ysize/2)
        elif croploc == 'upper right':
            xstart = xsizeorig - xsize + 1
            xstop = xsizeorig + 1
            ystart = ysizeorig - ysize
            ystop = ysizeorig + 1
        elif croploc == 'upper left':
            xstart = 1
            xstop = xsize + 1
            ystart = ysizeorig - ysize + 1
            ystop = ysizeorig + 1
        elif croploc == 'lower left':
            xstart = 1
            xstop = xsize + 1
            ystart = 1
            ystop = ysize + 1
        elif croploc == 'lower right':
            xstart = xsizeorig - xsize + 1
            xstop = xsizeorig + 1
            ystart = 1
            ystop = ysize + 1
        elif croploc == 'upper center':
            xcenter = int(xsizeorig/2)
            xstart = xcenter - int(xsize/2) + 1
            xstop = xcenter + int(xsize/2) + 1
            ystart = ysizeorig - ysize
            ystop = ysizeorig + 1
        elif croploc == 'centerlow':
            center = [int(xsizeorig/2), int(ysizeorig/2)]
            xstart = center[0] - int(xsize/2) + 1
            xstop = center[0] + int(xsize/2) + 1
            ystart = center[1]
            ystop = center[1] + int(ysize)
        else:
            print('croploc not recognized.')
            return None

        # Check that cropping dimensions are OK
        if any((xstart < 1, xstop < 1, ystart < 1,ystop < 1)):
            print('xsize/ysize dimensions are too large.')
            return None
        if any((xstart > xsizeorig+1, xstop > xsizeorig+1)):
            print('xsize dimensions are too large.')
            return None
        if any((ystart > ysizeorig+1, ystop > ysizeorig+1)):
            print('ysize dimensions are too large.')
            return None

        # Crop the image
        Im = Im[ystart:ystop, xstart-1:xstop]
        FitsHDU[0].data = Im

        # Write it to a new file
        if prefix is not None:
            tmpprefix = prefix
        else:
            tmpprefix = ''
        if suffix is not None:
            tmpsuffix = suffix
        else:
            tmpsuffix = ''

        OutFile = filepath + tmpprefix + filename + tmpsuffix + '.fits'
        if os.path.exists(OutFile):
            os.remove(OutFile)
        FitsHDU.writeto(OutFile)

    return None
