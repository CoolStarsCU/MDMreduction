#!/usr/bin/env python
#!/user/covey/iraf/mypython

###############################################################################
#spectra_splitter.py
#
# Takes out a specific order of a NIRSPEC IR spectra.  Call requires the
# name of the spectra the order will be extracted from, as well as the
# number of the order to extract from.
#
#Kevin Covey
#Version 1.1 (7-23-03)
###############################################################################

from pyraf import iraf
from pyraf.iraf import imutil

def spectra_splitter(image,order):

    orders = ['[1:1024,1:115]','[1:1024,116:260]','[1:1024,261:400]','[1:1024,401:545]','[1:1024,546:710]','[1:1024,711:850]','[1:1024,851:1024]']

    iraf.unlearn(iraf.imcopy)

    iraf.imcopy(input = image + orders[order-1] , output = 'o' + str(order) + '.' + image)
    
    iraf.flprcache()
    
