#!/usr/bin/env python
#!/user/covey/iraf/mypython

###############################################################################
# wavecal.py
#
# For MODspec spectra taken at MDM, this pyraf routine helps automate the wavelength
# calibration process, including the identification of lamp lines, re-identification
# of additional lamps, and application of wavelength shifts based on the 5577.339
# OH night sky line.
#
# This pipeline should only be run *after* running the MODprep.pro and raw2extract.py scripts.
#
# At that point, you should be ready to:
#
#   1. save a copy of the file you used to run raw2extract.py (e.g., data.tbl -> towavecal.tbl);
#      edit to only include lines listing lamps, science targets, and flux standards
#

#   11. After the spectra have been extracted, MODpipeline will use the first flux standard as a trace to extract the
#          lamps, and then ask you to identify all the lines.  use 'm' to mark each line with its wavelength,
#          then 'f' to fit the lines (tweak order with :o # and delete bad points with 'd').  'q' when done, both
#          with the fitting and the whole process.
#
#   12. MODpipeline will then apply that wavelength solution to the flux standards and science targets, linearize
#          the spectra, and use SETAIRMASS to get ready for flux calibration.

#Kevin Covey
#Version 1.0 (Last Modified 3-12-12; segment of former MODpipeline.py code)
#Stephanie Douglas
#Version 2.0 (Last Modified 4-23-15)
###############################################################################

import os

import numpy as np

from pyraf import iraf
from pyraf.iraf import noao
from pyraf.iraf import imutil, imred, crutil, ccdred, echelle, images, tv
from pyraf.iraf import system, twodspec, longslit, apextract, onedspec, astutil

from list_utils import read_reduction_list

def combine_lamps(hgne_lamp,xe_lamp,output_lamp):
    """ combine Xenon and HgNe lamps to use both sets of lines together"""
    print "combining HgNe: {} and Xe: {}".format(hgne_lamp,xe_lamp)

    xe_temp = "temp."+xe_lamp+".fits"
    # multiply the xenon lamp so that its lines show up stronger 
    # against the HgNe lines
    noao.onedspec.sarith(xe_lamp,"*",3.0,xe_temp)
    print "made new Xe lamp"
    noao.onedspec.sarith(xe_temp,"+",hgne_lamp,output_lamp)
    print "made new combination lamp"

    # now delete the temporary lamp
    os.remove(xe_temp)

    # return the name of the new lamp
    return output_lamp
    

def wavecal(imagelist,hgne_lamp=None,xe_lamp=None):

# get subsets of spectra: flats, lamps, biases, objects and std spectra

    image_dict = read_reduction_list(imagelist,obj_types=["obj","std","lamp"])

    # split out the relevant lists 
    # and make them arrays for marginally faster readout
    lamp_list = np.array(image_dict["lamp"])
    science_list = np.array(image_dict["science_list"])
    science_names = np.array(image_dict["science_names"])
    science_lamps = np.array(image_dict["science_reference_lamp"])
    std_list = np.array(image_dict["std"])
    std_names = np.array(image_dict["std_names"])

    numlamps = len(lamp_list)
    numscience = len(science_list)
    numstds = len(std_list)
    
# now identify lamp features
    
    # If an HgNe lamp and Xe lamp are provided, combine them
    if (hgne_lamp is not None) and (xe_lamp is not None):
        print "creating new reference lamp lamp.HgNeXe"
        new_ref_lamp = combine_lamps("lamp."+hgne_lamp,"lamp."+xe_lamp,
                                     "lamp.HgNeXe")
        reference_lamp = "HgNeXe"
    # If there's only an HgNe lamp, that's the reference lamp
    # (means that it's HgNe+Ne)
    elif (hgne_lamp is not None) and (xe_lamp is None):
        print "using HgNe reference lamp"
        reference_lamp = hgne_lamp
    # Otherwise just use the first lamp on reduction list
    else:
        reference_lamp = science_lamps[0]
        print "reference lamp is " + reference_lamp

    iraf.noao.onedspec.identify(images = 'lamp.'+reference_lamp)

# applying the lamp spectrum to the other lamps shouldn't matter, 
# because we're only using the reference lamp anyway

# now apply lamp to science objects

    iraf.noao.onedspec.refspec.sort = 'none'
    iraf.noao.onedspec.refspec.group = ''
    iraf.noao.onedspec.refspec.time = 'no'
    iraf.noao.onedspec.refspec.override = 'yes'
    iraf.noao.onedspec.refspec.confirm = 'no'
    iraf.noao.onedspec.refspec.assign = 'yes'
    iraf.noao.onedspec.dispcor.linearize = 'no'

    os.mkdir('wavecal')

    for j,science_file in enumerate(science_list):
        iraf.noao.onedspec.refspec(input = 'extract/ex.' + science_file, 
                                   reference = 'lamp.'+reference_lamp)
        iraf.noao.onedspec.dispcor(input = 'extract/ex.' + science_file, 
                                   output = 'wavecal/wc.' + science_file)
        iraf.noao.onedspec.scopy(input = 'wavecal/wc.' + science_file, 
                                 output = 'wavecal/wc.' + science_file, 
                                 bands = '3', format = "onedspec")
        iraf.noao.onedspec.wspectext(input = 'wavecal/wc.' + science_file + '.2001', 
                                     output = 'wavecal/sky.' + science_file,  
                                     header = 'no')

    iraf.flprcache()
