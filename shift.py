#!/usr/bin/env python
#!/user/covey/iraf/mypython

###############################################################################
# finalize.py
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
###############################################################################

from pyraf import iraf
from pyraf.iraf import noao
from pyraf.iraf import imutil, imred, crutil, ccdred, echelle, images, tv
from pyraf.iraf import system, twodspec, longslit, apextract, onedspec, astutil
import spectra_splitter, split_strings, sky_checker, os
from split_strings import split_strings
from spectra_splitter import spectra_splitter
from flat_normalizer import flat_normalizer
from sky_checker import sky_checker

import numpy as np

from list_utils import read_OI_shifts

def shift(imagelist):

# get subsets of spectra: flats, lamps, biases, objects and std spectra

    image_dict = read_OI_shifts(imagelist,science_types=["obj","std"])

    # split out the relevant lists 
    # and make them arrays for marginally faster readout
    science_list = np.array(image_dict["science_list"])
    science_names = np.array(image_dict["science_names"])
    shifts = np.array(image_dict["science_shift"])
    shift_errs = np.array(image_dict["science_shift_err"])
    shift_qual = np.array(image_dict["science_shift_qual"])
    std_list = np.array(image_dict["std"])
    std_names = np.array(image_dict["std_names"])

    numscience = len(science_list)
    numstds = len(std_list)


# now apply the shift determined by the IDL code to each spectrum, and add lines
# in the header to document the shift applied, the error in the shift, and the
# quality factor of the sky spectrum that was used to derive the shift..
    for j,science_file in enumerate(science_list):
        print 'wavecal/wc.'+science_file
        iraf.noao.onedspec.scopy(input = 'wavecal/wc.' + science_file, 
                                 output = 'wavecal/preshift.' + science_file,  
                                 bands = '1', format = "onedspec")
        iraf.noao.onedspec.wspectext(
                    input = 'wavecal/preshift.{}.2001'.format(science_file), 
                    output = 'wavecal/preshift.' + science_file, header = 'no')
        iraf.images.imutil.hedit(images = 'wavecal/wc.'+science_file,  
                                 fields = 'OI_SHIFT', value = shifts[j],
                                 add = 'yes', show = 'yes')
        iraf.images.imutil.hedit(images = 'wavecal/wc.'+science_file,   
                                 fields = 'OI_ERR', value = shift_errs[j],
                                 add = 'yes', show = 'yes')
        iraf.images.imutil.hedit(images = 'wavecal/wc.'+science_file,   
                                 fields = 'OI_QUAL', value = shift_qual[j], 
                                  add = 'yes', show = 'yes') 

        iraf.noao.onedspec.dispcor(input = 'wavecal/wc.'+science_file,  
                                   output = 'wavecal/dc.'+science_file,   
                                   linearize = 'yes')
        iraf.noao.onedspec.specshift(spectra = 'wavecal/dc.'+science_file,   
                                     shift = shifts[j])
        iraf.noao.onedspec.scopy(input = 'wavecal/wc.' + science_file,   
                                 output = 'wavecal/postshift.' + science_file, 
                                 bands = '1', format = "onedspec")
        iraf.noao.onedspec.wspectext(
                    input = 'wavecal/postshift.{}.2001'.format(science_file),
                    output = 'wavecal/postshift.' + science_file, header = 'no')


    iraf.flprcache()
