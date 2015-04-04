#!/usr/bin/env python
#!/user/covey/iraf/mypython

###############################################################################
#flat_normalizer.py
#
# read in a flat spectra, plot a histogram of all the flux values, define the
# boundary between order and interorder pixels, divide by the mean of the
# order pixels, and replace all pixels lower than the boundary / normalization
# with 1
#
#Kevin Covey
#Version 1.0 (7-20-03)
###############################################################################

from pyraf import iraf
from pyraf.iraf import images, imutil
import split_strings
from split_strings import split_strings
def flat_normalizer(inputimage):

# get rid of any parameters that are stuck in the packages I'll need

    iraf.unlearn(iraf.imhistogram)
    iraf.unlearn(iraf.imstat)
    iraf.unlearn(iraf.imarith)

# make a histogram of the pixel values in the flat

    iraf.imhistogram.hist_type = 'normal'
    iraf.imhistogram.listout = 'no'
    iraf.imhistogram.nbins = '100'
    iraf.imhistogram.z1 = '0'
    iraf.imhistogram.z2 = '30000'
    iraf.imhistogram(image = inputimage)

# showing the histogram (preferably as a plot, but possibly as text), ask user to pick
# boundary between interorder pixels and order pixels

    boundary = raw_input('Value seperating order and interorder pixels: ')

# find the mean of the order pixels

    flatstat = iraf.imstat(images = inputimage, fields = 'mean', lower = str(boundary), Stdout = 1)

    splitflat = split_strings(flatstat[1])

    mean = float(splitflat[0])

# divide the flat by the mean of order pixels
    
    outname = 'flats/flat.' + inputimage[0:]

    iraf.imarith(operand1 = inputimage, op = '/', operand2 = mean, result = outname)

# overwrite with a 1 interorder pixels

    getridof = float(boundary) / float(mean)

    iraf.imreplace(images = outname, upper = str(getridof), value = '1')
    
    iraf.flprcache()
