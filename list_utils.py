# utility for reading in lists of objects made by MODprep 
# (which I should convert to python...)

import astropy.io.ascii as at
import astropy.io.fits as fits
import numpy as np
import subprocess

def MODprep(ref_lamp, input_file="MODprep.lis", output_file="to_reduce.lis", 
            data_section="[1:400,1:2000]", bias_section="[401:464,1:2000]"):
    
    subprocess.call("ls *.fit > " + input_file, shell=True)
    
    infile = at.read(input_file, data_start=0)
    names = infile.columns[0].data
    
    with open(output_file,"w") as f:
        for i, name in enumerate(names):
            hdu = fits.open(name)
            imgtyp = hdu[0].header["IMAGETYP"].strip().upper()
            #datareg = hdu[0].header["DATASEC"].strip()
            #biasreg = hdu[0].header["BIASSEC"].strip()
            obj_name = hdu[0].header["OBJECT"].strip()
            
            if imgtyp in ["FOCUS", "LAMP", "LAMPS"]:
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
                    name[:-4], category, obj_name, data_section, 
                    bias_section, ref_lamp))
        f.close()

def read_list(imagelist,
              column_headers=["ccdno","type","target","image_region",
                              "bias_region","reference_lamp"],
              obj_types=["obj","std","flat","bias","lamp"],
              science_types=["obj","std"],return_other_cols=True,
              return_regions=True):
    """
    read in an object list made by MODprep. 

    input
    -----
    imagelist: string, filename
        filename for a list of spectra output by MODprep.pro

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
    image_list = at.read(imagelist,data_start=0,
                         names=column_headers)

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

def read_reduction_list(imagelist,obj_types=["obj","std","flat","bias","lamp"],
                        science_types=["obj","std"],return_regions=True):
    """ read the list of files to be reduced that was output by MODprep.pro """
    output = read_list(imagelist,obj_types=obj_types,science_types=science_types,
                       return_regions=return_regions,return_other_cols=True)
    junk = output.pop("science_image_region")
    junk = output.pop("science_bias_region")
    return output

def read_OI_shifts(shiftlist,science_types=["obj","std"]):
    """ read the list of files that was output by OIshift_corr.py
    """
    output = read_list(shiftlist,column_headers=["ccdno","type","target",
                                                 "shift","shift_err","shift_qual"],
                       obj_types=science_types,science_types=science_types,
                       return_regions=False,return_other_cols=True)
    return output

def generate_shift_list(imagelist,output_filename,science_types=["obj","std"]):
    """ 
    takes the original list created by MODprep() and removes all calibration
    files. The output list can be used to run OIshift_corr.pro and fluxcal.py
    """

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
    """ 
    make sure all names in the output target list are unique 
    input
    -----
    science_names: array_like containing strings

    output
    ------
    duplicates: bool
        whether the input list contains duplicated names
    """

    duplicates = True
    
    unique_names,name_counts = np.unique(science_names,return_counts=True)

    if len(unique_names)==len(science_names):
        print("All target names are unique!")
        duplicates = False
    else:
        print("The following names are repeated:")
        print(unique_names[name_counts>1])

    return duplicates
