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

def wavecal(imagelist):

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
    lamplist = []
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
            elif spectra[j][1] == 'lamp':
                lamplist.append(spectra[j])
                j = j + 1
            else:
                j = j + 1

# find the number of objects in each list

    numlamps = len(lamplist)
    numscience = len(sciencelist)
    numstds = len(stdlist)
    
# now identify lamp features

    iraf.noao.onedspec.identify(images = 'lamp.'+lamplist[0][0])

# now apply lamp solution to the rest of the lamps

    k = 0
    while k < numlamps:
         iraf.noao.onedspec.reidentify(reference = 'lamp.'+lamplist[0][0], images = 'lamp.'+lamplist[k][0])
         k = k+1

# now apply the appropriate lamp to each spectrum

# now apply lamp to science objects

    iraf.noao.onedspec.refspec.sort = 'none'
    iraf.noao.onedspec.refspec.group = ''
    iraf.noao.onedspec.refspec.time = 'no'
    iraf.noao.onedspec.refspec.override = 'yes'
    iraf.noao.onedspec.refspec.confirm = 'no'
    iraf.noao.onedspec.refspec.assign = 'yes'
    iraf.noao.onedspec.dispcor.linearize = 'no'

    os.mkdir('wavecal')

    e = 0
    while e < numscience:
        #print sciencelist[e][5]
        iraf.noao.onedspec.refspec(input = 'extract/ex.' + sciencelist[e][0], reference = 'lamp.'+sciencelist[e][5])
        iraf.noao.onedspec.dispcor(input = 'extract/ex.' + sciencelist[e][0], output = 'wavecal/wc.' + sciencelist[e][0])
        iraf.noao.onedspec.scopy(input = 'wavecal/wc.' + sciencelist[e][0], output = 'wavecal/wc.' + sciencelist[e][0], bands = '3', format = "onedspec")
        iraf.noao.onedspec.wspectext(input = 'wavecal/wc.' + sciencelist[e][0] + '.2001', output = 'wavecal/sky.' + sciencelist[e][0], header = 'no')
        e = e + 1 

#    f = 0
#    while f < numscience:


    iraf.flprcache()
