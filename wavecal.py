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

def wavecal(imagelist):

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

    reference_lamp = lamp_list[0]

    iraf.noao.onedspec.identify(images = 'lamp.'+reference_lamp)

# now apply lamp solution to the rest of the lamps
# not that it matters since we apply the first lamp to everything anyway

    for lamp_file in lamp_list:
         iraf.noao.onedspec.reidentify(reference = 'lamp.'+reference_lamp,  
                                       images = 'lamp.'+lamp_file)

# now apply the appropriate lamp to each spectrum - missing?

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
                                   reference = 'lamp.'+science_lamps[j])
        iraf.noao.onedspec.dispcor(input = 'extract/ex.' + science_file, 
                                   output = 'wavecal/wc.' + science_file)
        iraf.noao.onedspec.scopy(input = 'wavecal/wc.' + science_file, 
                                 output = 'wavecal/wc.' + science_file, 
                                 bands = '3', format = "onedspec")
        iraf.noao.onedspec.wspectext(input = 'wavecal/wc.' + science_file + '.2001', 
                                     output = 'wavecal/sky.' + science_file,  
                                     header = 'no')

    iraf.flprcache()
