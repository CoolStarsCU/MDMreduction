#!/usr/bin/env python
#!/user/covey/iraf/mypython

import os
import numpy as np

from pyraf import iraf
from pyraf.iraf import noao
from pyraf.iraf import imutil, imred, crutil, ccdred, echelle, images, tv
from pyraf.iraf import system, twodspec, longslit, apextract, onedspec, astutil

from list_utils import read_OI_shifts

def shift(imagelist="OI_shifts.tbl"):

    '''
    This function applies wavelength shifts based on the 5577.339 OH night sky line.

    imagelist - String; filename of list of spectra to process (output from OI_shiftcorr.main())

    Kevin Covey
    version: 1.0 (Last Modified 3-12-12; segment of former MODpipeline.py code)
    Stephanie Douglas
    version 2.0 (Last Modified 4-23-15)
    Alejandro Nunez
    version 2.1 (Last Modified 2017-11_)
    '''
    # Read input file
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


# now apply the shift determined by the code to each spectrum, and add lines
# in the header to document the shift applied, the error in the shift, and the
# quality factor of the sky spectrum that was used to derive the shift..
    for j,science_file in enumerate(science_list):
        print 'wavecal/wc.'+science_file
        iraf.noao.onedspec.scopy(input = 'wavecal/wc.' + science_file,
                                 output = 'wavecal/preshift.' + science_file,
                                 bands = '1', format = "onedspec")
        iraf.noao.onedspec.wspectext(
                    input = 'wavecal/preshift.{}.0001'.format(science_file),
                    output = 'wavecal/preshift.' + science_file, header = 'no')
        iraf.images.imutil.hedit(images = 'wavecal/wc.'+science_file,
                                 fields = 'OI_SHIFT', value = shifts[j],
                                 add = 'yes', show = 'yes', verify='no')
        iraf.images.imutil.hedit(images = 'wavecal/wc.'+science_file,
                                 fields = 'OI_ERR', value = shift_errs[j],
                                 add = 'yes', show = 'yes', verify='no')
        iraf.images.imutil.hedit(images = 'wavecal/wc.'+science_file,
                                 fields = 'OI_QUAL', value = shift_qual[j],
                                  add = 'yes', show = 'yes', verify='no')
        iraf.noao.onedspec.dispcor(input = 'wavecal/wc.'+science_file,
                                   output = 'wavecal/dc.'+science_file,
                                   linearize = 'yes')
        iraf.noao.onedspec.specshift(spectra = 'wavecal/dc.'+science_file,
                                     shift = shifts[j])
        iraf.noao.onedspec.scopy(input = 'wavecal/wc.' + science_file,
                                 output = 'wavecal/postshift.' + science_file,
                                 bands = '1', format = "onedspec")
        # Obsolete: (?)
        # iraf.noao.onedspec.wspectext(
        #             input = 'wavecal/postshift.{}.2001'.format(science_file),
        #             output = 'wavecal/postshift.' + science_file, header = 'no')


    iraf.flprcache()
