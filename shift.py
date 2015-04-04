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

def shift(imagelist):

# read in a list of spectra as defined in point 1 above:

    input = open(imagelist, 'r')
    spectra = input.readlines()

# delete first two lines of the table, assuming a line of column headers
# and a line of whitespace
   # del spectra[0]
   # del spectra[0]

# count number of spectra to reduce
    spectranum = len(spectra)

# go through table line by line and remove spaces, tabs, and end of line marks
# each cleaned line is turned into a list, where spaces and tabs separate
# items in a list for each line.  Each line's list is itself nested in a larger
# list

    i = 0
    while i < spectranum:
        spectra[i] = split_strings(spectra[i])
        i = i + 1
        
# now create important subsets: flats, lamps, biases, objects and std spectra

#first make lists for all the objects
    sciencelist = []
    stdlist = []
    
#fill those lists with the names of objects

    j = 0
    while j < spectranum:
            if spectra[j][1] == 'obj':
                sciencelist.append(spectra[j])
                j = j + 1
            elif spectra[j][1] == 'std':
                sciencelist.append(spectra[j])
                stdlist.append(spectra[j])
                j = j + 1
            else:
                j = j + 1

# find the number of objects in each list

    numscience = len(sciencelist)
    numstds = len(stdlist)
    
# now apply the shift determined by the IDL code to each spectrum, and add lines
# in the header to document the shift applied, the error in the shift, and the
# quality factor of the sky spectrum that was used to derive the shift..
    e = 0
    while e < numscience:

        print 'wavecal/wc.'+sciencelist[e][0]
        iraf.noao.onedspec.scopy(input = 'wavecal/wc.' + sciencelist[e][0], output = 'wavecal/preshift.' + sciencelist[e][0], bands = '1', format = "onedspec")
        iraf.noao.onedspec.wspectext(input = 'wavecal/preshift.' + sciencelist[e][0] + '.2001', output = 'wavecal/preshift.' + sciencelist[e][0], header = 'no')
        iraf.images.imutil.hedit(images = 'wavecal/wc.'+sciencelist[e][0], fields = 'OI_SHIFT', value = sciencelist[e][3], add = 'yes', show = 'yes')
        iraf.images.imutil.hedit(images = 'wavecal/wc.'+sciencelist[e][0], fields = 'OI_ERR', value = sciencelist[e][4], add = 'yes', show = 'yes')
        iraf.images.imutil.hedit(images = 'wavecal/wc.'+sciencelist[e][0], fields = 'OI_QUAL', value = sciencelist[e][5], add = 'yes', show = 'yes')                            
        iraf.noao.onedspec.dispcor(input = 'wavecal/wc.'+sciencelist[e][0], output = 'wavecal/dc.'+sciencelist[e][0], linearize = 'yes')
        iraf.noao.onedspec.specshift(spectra = 'wavecal/dc.'+sciencelist[e][0], shift = sciencelist[e][3])
        iraf.noao.onedspec.scopy(input = 'wavecal/wc.' + sciencelist[e][0], output = 'wavecal/postshift.' + sciencelist[e][0], bands = '1', format = "onedspec")
        iraf.noao.onedspec.wspectext(input = 'wavecal/postshift.' + sciencelist[e][0] + '.2001', output = 'wavecal/postshift.' + sciencelist[e][0], header = 'no')
        e = e + 1
    
    iraf.flprcache()
