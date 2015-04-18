#!/usr/bin/env python
#!/user/covey/iraf/mypython

###############################################################################
# raw2extract.py
#
# For MODspec spectra taken at MDM, this pyraf routine helps automate the reduction
# process from raw images through extracted spectra
#
# A typical useage case for this pipeline is:
#
#   1. Generate a file listing all your raw spectra.  In the directory with all your observations, run:
#                       ls *fits > raw_spectra.lis
#
#   2. Use the IDL program MODprep.pro to pull several useful parameters out of each stars header.
#         In IDL, simply run MODprep.pro; when prompted, specify the name of the file listing the
#         raw spectra, and the file MODprep should write its output to, such as prep_spectra.lis 
#         The output file will have 5 columns: for each spectra, the columns will contain the
#         filename (1), observation type [std = flux std, obj = science object, flat = flat,
#         bias = bias, wavelength lamps = lamp] (2), object name (3), data section of the ccd (4),
#         and overscan region (5).
#
#   3. Edit MODprep's output file to:
#         * find your flux standards, and change their observation type (in the second column) from 'obj' to 'std'
#         * if an object was observed multiple times, it will (initially) produce multiple spectra: apply
#           suffixes to the object names in the third column to distinguish them from one another
#
#   4. Confirm that you have the raw2extract.py file, as well as those it depends on (sky_checker.py, split_strings.py,
#         spectra_splitter.py, flat_normalizer.py) in your path or in the directory with your data
#
#   5. start pyraf in an xgterm window (make sure pyraf has access to your iraf login.cl file too)
#
#   6. import and call the reduction pipeline with:
#           import raw2extract
#           raw2extract.raw2extract('prep_spectra.lis') <--- replacing prep_spectra.lis with your file name
#
#   7. raw2extract will happily run in the background for a bit; it will keep busy by:
#         * combining the biases into masterbias.fits
#         * run ccd proc on all images to apply overscan correction, master bias
#           correction, and trim images
#         * clean science images with cosmicrays
#         * combine flats into masterflat.fits
#
#   8. The user will then be asked to help raw2extract to run the RESPONSE task to remove
#          the large-scale structure in the master flat.  The user should verify that the fit describes
#          the structure in the flat well; if the fit is poor, the order can be increased (or decreased) with
#          a :o # command (with the # giving the new order of the fit), followed by an 'f' and 'r' command to
#          redo the fit and redraw the results.
#
#   9. Once the masterflat has been normalized by the response, raw2extract will use normflat.fits to
#          flat-divide each (2-D) spectrum with ccdproc.
#
#   10. raw2extract will then ask the user to help extract the science targets and flux standards.
#          This will be done with standard apall tasks -- for more on this, see "A User's Guide to
#          Reducing Slit Spectra with IRAF" by Phil Massey.
#
#          ** Quick hints:
#
#            * To tweak background, hit b to enter the
#              background fitting window.  Then you can delete bad background sampling regions with 'z',
#              outline new regions with 's', redo the fit with 'f', and redraw the graph with 'r'.
#
#            * For some reason I can't figure out, some interactive y/n questions can only be answered in the
#              graphics window or the terminal window.  This means you have to mouse between them once for each
#              spectrum.  Annoying, but I can't find another solution. 
#
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

def raw2extract(imagelist):

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
        #print i, spectranum
        spectra[i] = split_strings(spectra[i])
        i = i + 1
        
# now create important subsets: flats, lamps, biases, objects and std spectra

#first make lists for all the objects
    flatlist = [] 
    lamplist = []
    biaslist = []
    sciencelist = []
    stdlist = []
    
#fill those lists with the names of objects

    j = 0
    while j < spectranum:
            #print j,len(spectra[j])
            if spectra[j][1] == 'obj':
                sciencelist.append(spectra[j])
                j = j + 1
            elif spectra[j][1] == 'std':
                sciencelist.append(spectra[j])
                stdlist.append(spectra[j])
                j = j + 1
            elif spectra[j][1] == 'flat':
                flatlist.append(spectra[j])
                j = j + 1
            elif spectra[j][1] == 'lamp':
                lamplist.append(spectra[j])
                j = j + 1
            elif spectra[j][1] == 'bias':
                biaslist.append(spectra[j])
                j = j + 1
            else:
                j = j + 1

# find the number of objects in each list

    numflats = len(flatlist)
    numlamps = len(lamplist)
    numscience = len(sciencelist)
    numbiases = len(biaslist)
    numstds = len(stdlist)
    
    overscanregion = flatlist[0][4] #'[' + redlowcoloverscan + ':' + redhighcoloverscan + ',' + redlowlinedata + ':' + redhighlinedata + ']'
    gooddata = flatlist[0][3]#'[' + redlowcoldata + ':' + redhighcoldata + ',' + redlowlinedata + ':' + redhighlinedata +']'

#define the gain and read noise
    gain = '2.7' # taken from header in Dec. 2010, but matches 1997 docs.
    readnoise = '7.9' # taken from web documentation dated 1997

#combine biases into master biases of each color

    k = 0
    biasestocombine = ''
    while k < numbiases:
        if k == 0:
            biasestocombine = biaslist[k][0]
            k = k + 1
        else:
            biasestocombine = biasestocombine + ',' + biaslist[k][0]
            k = k + 1

    iraf.noao.imred.ccdred.zerocombine.combine = 'median'
    iraf.noao.imred.ccdred.zerocombine.reject = 'minmax'
    iraf.noao.imred.ccdred.zerocombine.ccdtype = ''
    iraf.noao.imred.ccdred.zerocombine.process = 'no'
    iraf.noao.imred.ccdred.zerocombine.delete = 'no'
    iraf.noao.imred.ccdred.zerocombine.clobber = 'no'
    iraf.noao.imred.ccdred.zerocombine.scale = 'no'
    iraf.noao.imred.ccdred.zerocombine.statsec = ''
    iraf.noao.imred.ccdred.zerocombine.nlow = '0'
    iraf.noao.imred.ccdred.zerocombine.nhigh = '1'
    iraf.noao.imred.ccdred.zerocombine.nkeep = '1'
    iraf.noao.imred.ccdred.zerocombine.mclip = 'yes'
    iraf.noao.imred.ccdred.zerocombine.lsigma = '3.0'
    iraf.noao.imred.ccdred.zerocombine.hsigma = '3.0'
    iraf.noao.imred.ccdred.zerocombine.rdnoise = readnoise
    iraf.noao.imred.ccdred.zerocombine.gain = gain
    iraf.noao.imred.ccdred.zerocombine.snoise = '0.' #last thing that should change
    iraf.noao.imred.ccdred.zerocombine.pclip = '-0.5'
    iraf.noao.imred.ccdred.zerocombine.blank = '0.0'

    iraf.noao.imred.ccdred.zerocombine(input = biasestocombine, output = 'masterbias')

#now actually trim all the images and store them in a trimmed, overscan region

    iraf.noao.imred.ccdred.ccdproc.ccdtype = ''
    iraf.noao.imred.ccdred.ccdproc.noproc = 'no'
    iraf.noao.imred.ccdred.ccdproc.fixpix = 'no'
    iraf.noao.imred.ccdred.ccdproc.overscan = 'yes'
    iraf.noao.imred.ccdred.ccdproc.trim = 'yes'
    iraf.noao.imred.ccdred.ccdproc.zerocor = 'yes'
    iraf.noao.imred.ccdred.ccdproc.darkcor = 'no'
    iraf.noao.imred.ccdred.ccdproc.flatcor = 'no'
    iraf.noao.imred.ccdred.ccdproc.illumcor = 'no'
    iraf.noao.imred.ccdred.ccdproc.fringecor = 'no'
    iraf.noao.imred.ccdred.ccdproc.readcor = 'no'
    iraf.noao.imred.ccdred.ccdproc.scancor = 'no'
    iraf.noao.imred.ccdred.ccdproc.readaxis = 'line'
    iraf.noao.imred.ccdred.ccdproc.interactive = 'no'
    iraf.noao.imred.ccdred.ccdproc.function = 'legendre'
    iraf.noao.imred.ccdred.ccdproc.order = '1'
    iraf.noao.imred.ccdred.ccdproc.sample = '*'
    iraf.noao.imred.ccdred.ccdproc.naverage = '1'
    iraf.noao.imred.ccdred.ccdproc.niterate = '1'
    iraf.noao.imred.ccdred.ccdproc.low_reject = '3.0'
    iraf.noao.imred.ccdred.ccdproc.high_reject = '3.0'

    os.mkdir('trimmed')

    iraf.noao.imred.ccdred.ccdproc.biassec = overscanregion
    iraf.noao.imred.ccdred.ccdproc.trimsec = gooddata
    iraf.noao.imred.ccdred.ccdproc.zero = 'masterbias'
   
    m = 0
    while m < numscience:
        iraf.noao.imred.ccdred.ccdproc(images = sciencelist[m][0], output = 'trimmed/tr.' + sciencelist[m][0])
        m = m + 1

    n = 0
    while n < numflats:
        iraf.noao.imred.ccdred.ccdproc(images = flatlist[n][0], output = 'trimmed/tr.' + flatlist[n][0])
        n = n + 1

    o = 0
    while o < numlamps:
        iraf.noao.imred.ccdred.ccdproc(images = lamplist[o][0], output = 'trimmed/tr.' + lamplist[o][0])
        o = o + 1

    p = 0
    while p < numbiases:
        iraf.noao.imred.ccdred.ccdproc(images = biaslist[p][0], output = 'trimmed/tr.' + biaslist[p][0])
        p = p + 1    

# now clean out the cosmic rays (just from the science images)

    os.mkdir('cleaned')

    iraf.unlearn(iraf.imred.crutil.cosmicrays)

    iraf.imred.crutil.cosmicrays.threshold = '25'
    iraf.imred.crutil.cosmicrays.fluxratio = '2'
    iraf.imred.crutil.cosmicrays.npasses = '5'
    iraf.imred.crutil.cosmicrays.window = '5'
    iraf.imred.crutil.cosmicrays.interactive = 'no'

    print numscience
    u = 0
    while u < numscience:
        iraf.imred.crutil.cosmicrays(input = 'trimmed/tr.' + sciencelist[u][0], output = 'cleaned/cr.' + sciencelist[u][0])
        u = u + 1
    
# now combine the flats and use response on the master flats

    w = 0
    flatstocombine = ''
    while w < numflats:
        if w == 0:
            flatstocombine = 'trimmed/tr.' + flatlist[w][0]
            w = w + 1
        else:
            flatsstocombine = flatstocombine + ',' + 'trimmed/tr.' + flatlist[w][0]
            w = w + 1

    iraf.noao.imred.ccdred.flatcombine.combine = 'average'
    iraf.noao.imred.ccdred.flatcombine.reject = 'avsigclip'
    iraf.noao.imred.ccdred.flatcombine.ccdtype = ''
    iraf.noao.imred.ccdred.flatcombine.process = 'no'
    iraf.noao.imred.ccdred.flatcombine.subsets = 'no'
    iraf.noao.imred.ccdred.flatcombine.delete = 'no'
    iraf.noao.imred.ccdred.flatcombine.clobber = 'no'
    iraf.noao.imred.ccdred.flatcombine.scale = 'mode'
    iraf.noao.imred.ccdred.flatcombine.nlow = '0'
    iraf.noao.imred.ccdred.flatcombine.nhigh = '1'
    iraf.noao.imred.ccdred.flatcombine.nkeep = '1'
    iraf.noao.imred.ccdred.flatcombine.mclip = 'yes'
    iraf.noao.imred.ccdred.flatcombine.lsigma = '3.0'
    iraf.noao.imred.ccdred.flatcombine.hsigma = '3.0'
    iraf.noao.imred.ccdred.flatcombine.rdnoise = readnoise
    iraf.noao.imred.ccdred.flatcombine.gain = gain
    iraf.noao.imred.ccdred.flatcombine.snoise = '0.' #last thing that should change
    iraf.noao.imred.ccdred.flatcombine.pclip = '-0.5'
    iraf.noao.imred.ccdred.flatcombine.blank = '0.0'

    iraf.noao.imred.ccdred.flatcombine(input = flatstocombine, output = 'masterflat')

    iraf.noao.twodspec.longslit.dispaxis = '2'
    iraf.noao.twodspec.longslit.response.interactive = 'yes'
    iraf.noao.twodspec.longslit.response.threshold = 'INDEF'
    iraf.noao.twodspec.longslit.response.sample = '*'
    iraf.noao.twodspec.longslit.response.naverage = '1'
    iraf.noao.twodspec.longslit.response.function = 'spline3'
    iraf.noao.twodspec.longslit.response.order = '11'
    iraf.noao.twodspec.longslit.response.low_reject = '3.'
    iraf.noao.twodspec.longslit.response.high_reject = '3.'
    iraf.noao.twodspec.longslit.response.niterate = '1'
    iraf.noao.twodspec.longslit.response.grow = '0'

    iraf.noao.twodspec.longslit.response(calibration = 'masterflat', normalization = 'masterflat', response = 'normflat') 

# now actually do the flat division

    iraf.noao.imred.ccdred.ccdproc.ccdtype = ''
    iraf.noao.imred.ccdred.ccdproc.noproc = 'no'
    iraf.noao.imred.ccdred.ccdproc.fixpix = 'no'
    iraf.noao.imred.ccdred.ccdproc.overscan = 'no'
    iraf.noao.imred.ccdred.ccdproc.trim = 'no'
    iraf.noao.imred.ccdred.ccdproc.zerocor = 'no'
    iraf.noao.imred.ccdred.ccdproc.darkcor = 'no'
    iraf.noao.imred.ccdred.ccdproc.flatcor = 'yes'
    iraf.noao.imred.ccdred.ccdproc.illumcor = 'no'
    iraf.noao.imred.ccdred.ccdproc.fringecor = 'no'
    iraf.noao.imred.ccdred.ccdproc.readcor = 'no'
    iraf.noao.imred.ccdred.ccdproc.scancor = 'no'
    iraf.noao.imred.ccdred.ccdproc.readaxis = 'line'
    iraf.noao.imred.ccdred.ccdproc.interactive = 'no'
    iraf.noao.imred.ccdred.ccdproc.function = 'legendre'
    iraf.noao.imred.ccdred.ccdproc.order = '1'
    iraf.noao.imred.ccdred.ccdproc.sample = '*'
    iraf.noao.imred.ccdred.ccdproc.naverage = '1'
    iraf.noao.imred.ccdred.ccdproc.niterate = '1'
    iraf.noao.imred.ccdred.ccdproc.low_reject = '3.0'
    iraf.noao.imred.ccdred.ccdproc.high_reject = '3.0'

    os.mkdir('flattened')

    iraf.noao.imred.ccdred.ccdproc.flat = 'normflat'
    
    y = 0
    while y < numscience:
        iraf.noao.imred.ccdred.ccdproc(images = 'cleaned/cr.' + sciencelist[y][0], output = 'flattened/fl.' + sciencelist[y][0])
        #iraf.imcopy('cleaned/cr.' + sciencelist[y][0],'flattened/fl.' + sciencelist[y][0])
        y = y + 1

    iraf.unlearn(iraf.noao.twodspec.apextract.apall)
    
    iraf.noao.twodspec.apextract.apall.apertures = '1'
    iraf.noao.twodspec.apextract.apall.format = 'multispec'
    iraf.noao.twodspec.apextract.apall.references = ''
    iraf.noao.twodspec.apextract.apall.profiles = ''
    iraf.noao.twodspec.apextract.apall.interactive = 'yes'#'no'#
    iraf.noao.twodspec.apextract.apall.find = 'yes'
    iraf.noao.twodspec.apextract.apall.recenter = 'yes'
    iraf.noao.twodspec.apextract.apall.resize = 'yes'
    iraf.noao.twodspec.apextract.apall.edit = 'yes'
    iraf.noao.twodspec.apextract.apall.trace = 'yes'
    iraf.noao.twodspec.apextract.apall.fittrace = 'yes'
    iraf.noao.twodspec.apextract.apall.extract = 'yes'
    iraf.noao.twodspec.apextract.apall.extras = 'yes'
    iraf.noao.twodspec.apextract.apall.review = 'no'
    iraf.noao.twodspec.apextract.apall.line = 'INDEF'
    iraf.noao.twodspec.apextract.apall.nsum = '10'
    iraf.noao.twodspec.apextract.apall.lower = '-5'
    iraf.noao.twodspec.apextract.apall.upper = '5'
    iraf.noao.twodspec.apextract.apall.b_function = 'chebyshev'
    iraf.noao.twodspec.apextract.apall.b_order = '2'
    iraf.noao.twodspec.apextract.apall.b_sample = '-30:-18,18:30'
    iraf.noao.twodspec.apextract.apall.b_naverage = '-100'
    iraf.noao.twodspec.apextract.apall.b_niterate = '3'
    iraf.noao.twodspec.apextract.apall.b_low_reject = '3'
    iraf.noao.twodspec.apextract.apall.b_high_reject = '3'
    iraf.noao.twodspec.apextract.apall.b_grow =  '0'
    iraf.noao.twodspec.apextract.apall.width = '10'
    iraf.noao.twodspec.apextract.apall.radius = '10'
    iraf.noao.twodspec.apextract.apall.threshold = '0'
    iraf.noao.twodspec.apextract.apall.nfind = '1'
    iraf.noao.twodspec.apextract.apall.minsep = '5'
    iraf.noao.twodspec.apextract.apall.maxsep = '1000'
    iraf.noao.twodspec.apextract.apall.order = 'increasing'
    iraf.noao.twodspec.apextract.apall.aprecenter = ''
    iraf.noao.twodspec.apextract.apall.npeaks = 'INDEF'
    iraf.noao.twodspec.apextract.apall.shift = 'no'
    iraf.noao.twodspec.apextract.apall.llimit = 'INDEF'
    iraf.noao.twodspec.apextract.apall.ulimit = 'INDEF'
    iraf.noao.twodspec.apextract.apall.ylevel = '.2'
    iraf.noao.twodspec.apextract.apall.peak = 'yes'
    iraf.noao.twodspec.apextract.apall.bkg = 'yes'
    iraf.noao.twodspec.apextract.apall.r_grow = '0.'
    iraf.noao.twodspec.apextract.apall.avglimits =  'no'  
    iraf.noao.twodspec.apextract.apall.t_nsum = '10'
    iraf.noao.twodspec.apextract.apall.t_step = '10'
    iraf.noao.twodspec.apextract.apall.t_nlost = '5'
    iraf.noao.twodspec.apextract.apall.t_function = 'legendre'
    iraf.noao.twodspec.apextract.apall.t_order = '5'
    iraf.noao.twodspec.apextract.apall.t_sample = '*'
    iraf.noao.twodspec.apextract.apall.t_naverage = '1'
    iraf.noao.twodspec.apextract.apall.t_niterate = '3'
    iraf.noao.twodspec.apextract.apall.t_low_reject = '3'
    iraf.noao.twodspec.apextract.apall.t_high_reject = '3'
    iraf.noao.twodspec.apextract.apall.t_grow = '0'
    iraf.noao.twodspec.apextract.apall.background = 'fit'
    iraf.noao.twodspec.apextract.apall.skybox = '1'
    iraf.noao.twodspec.apextract.apall.weights = 'variance'
    iraf.noao.twodspec.apextract.apall.pfit = 'fit1d'
    iraf.noao.twodspec.apextract.apall.clean = 'yes'
    iraf.noao.twodspec.apextract.apall.saturation = 'INDEF' #needs to be confirmed
    iraf.noao.twodspec.apextract.apall.readnoise = readnoise # taken from web documentation dated 1997
    iraf.noao.twodspec.apextract.apall.gain = gain # taken from header in Dec. 2010, but matches 1997 docs.
    iraf.noao.twodspec.apextract.apall.lsigma = 4
    iraf.noao.twodspec.apextract.apall.usigma = 4
    iraf.noao.twodspec.apextract.apall.nsubaps = '1'

    os.mkdir('extract')

    a = 0
    while a < numscience:
        iraf.noao.twodspec.apextract.apall(input = 'flattened/fl.' + sciencelist[a][0], output = 'extract/ex.' + sciencelist[a][0])
        a = a + 1

# now extract a blue and red lamp using the first science object as our trace.

    iraf.noao.twodspec.apextract.apall.references = 'flattened/fl.' + stdlist[0][0]
    iraf.noao.twodspec.apextract.apall.interactive = 'no'
    iraf.noao.twodspec.apextract.apall.find = 'no'
    iraf.noao.twodspec.apextract.apall.recenter = 'no'
    iraf.noao.twodspec.apextract.apall.resize = 'no'
    iraf.noao.twodspec.apextract.apall.edit = 'no'
    iraf.noao.twodspec.apextract.apall.trace = 'no'
    iraf.noao.twodspec.apextract.apall.fittrace = 'no'
    iraf.noao.twodspec.apextract.apall.extract = 'yes'
    iraf.noao.twodspec.apextract.apall.background = 'none'
    iraf.noao.twodspec.apextract.apall.clean = 'no'
    iraf.noao.twodspec.apextract.apall.weights = 'none'

    a = 0
    while a < numlamps: 
        iraf.noao.twodspec.apextract.apall(input = 'trimmed/tr.' + lamplist[a][0], output = 'lamp.'+lamplist[a][0])
        a = a + 1

    iraf.flprcache()
