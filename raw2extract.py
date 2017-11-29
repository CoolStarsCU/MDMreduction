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
import pdb

def raw2extract(instrument, imagelist='to_reduce.lis', apall_interactive='yes'):
    '''
    Instrument can be Modspec, OSMOS4k, or OSMOSr4k.
    '''
    # Get subsets of spectra: flats, lamps, biases, objects and std spectra
    # If instrument is OSMOS4k, then the bias correction steps should be done in
    # OSMOS_bias.py before running this.
    # If instrument is OSMOS r4k, then the bias correction steps should be done
    # using the script proc4k.py before running this.

    image_dict = read_reduction_list(imagelist)

    # Split out the relevant lists
    # and make them arrays for marginally faster readout
    flat_list = np.array(image_dict["flat"])
    lamp_list = np.array(image_dict["lamp"])
    if instrument.upper() == 'MODSPEC':
        bias_list = np.array(image_dict["bias"])
    science_list = np.array(image_dict["science_list"])
    science_names = np.array(image_dict["science_names"])
    std_list = np.array(image_dict["std"])
    std_names = np.array(image_dict["std_names"])

    # Find the number of objects in each list
    numflats = len(flat_list)
    numlamps = len(lamp_list)
    numscience = len(science_list)
    if instrument.upper() == 'MODSPEC':
        numbiases = len(bias_list)
    numstds = len(std_list)

    overscanregion = image_dict["overscan_region"]
    gooddata = image_dict["good_region"]

    print "{} Science images".format(numscience)

    # Define gain and read noise
    if instrument.upper() == 'MODSPEC':
        gain = '2.7' # taken from header in Dec. 2010, but matches 1997 docs.
    elif instrument.upper() == 'OSMOS4K':
        gain = '2.7'
    elif instrument.upper() == 'OSMOSR4K':
        gain = '1.0' # Unknown as of Nov 2017

    readnoise = '7.9' # taken from web documentation dated 1997

    # Combine biases into master biases of each color
    if instrument.upper() == 'MODSPEC':
        biasestocombine = ''
        for bias_file in bias_list:
            biasestocombine = biasestocombine + ", {}".format(bias_file)
        # Remove that first comma
        biasestocombine = biasestocombine[2:]

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

    # Trim all the images and store them in a trimmed, overscan region

    iraf.noao.imred.ccdred.ccdproc.ccdtype = ''
    iraf.noao.imred.ccdred.ccdproc.noproc = 'no'
    iraf.noao.imred.ccdred.ccdproc.fixpix = 'no'
    if 'OSMOS' in instrument.upper():
        iraf.noao.imred.ccdred.ccdproc.overscan = 'no'
        iraf.noao.imred.ccdred.ccdproc.zerocor = 'no'
    else:
        iraf.noao.imred.ccdred.ccdproc.overscan = 'yes'
        iraf.noao.imred.ccdred.ccdproc.zerocor = 'yes'
    iraf.noao.imred.ccdred.ccdproc.trim = 'yes'
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

    for science_file in science_list:
        iraf.noao.imred.ccdred.ccdproc(images = science_file,
                                       output = 'trimmed/tr.' + science_file)

    for flat_file in flat_list:
        iraf.noao.imred.ccdred.ccdproc(images = flat_file,
                                       output = 'trimmed/tr.' + flat_file)

    for lamp_file in lamp_list:
        iraf.noao.imred.ccdred.ccdproc(images = lamp_file,
                                       output = 'trimmed/tr.' + lamp_file)

    # Clean out the cosmic rays (just from the science images)
    os.mkdir('cleaned')

    iraf.unlearn(iraf.imred.crutil.cosmicrays)

    iraf.imred.crutil.cosmicrays.threshold = '25'
    iraf.imred.crutil.cosmicrays.fluxratio = '2'
    iraf.imred.crutil.cosmicrays.npasses = '5'
    iraf.imred.crutil.cosmicrays.window = '5'
    iraf.imred.crutil.cosmicrays.interactive = 'no'

    for science_file in science_list:
        print "cleaning " + science_file
        iraf.imred.crutil.cosmicrays(input = 'trimmed/tr.' + science_file,
                                     output = 'cleaned/cr.' + science_file)

    # Combine the flats and use response on the master flats

    flatstocombine = ''
    for flat_file in flat_list:
        flatstocombine = flatstocombine + ", {}".format(flat_file)
    # Remove that first comma
    flatstocombine = flatstocombine[2:]

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

    if 'OSMOS' in instrument.upper():
        iraf.noao.twodspec.longslit.dispaxis = '1' # horizontal dispersion
    else:
        iraf.noao.twodspec.longslit.dispaxis = '2' # vertical dispersion
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

    # Do the flat division
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
    if 'OSMOS' in instrument.upper():
        iraf.noao.imred.ccdred.ccdproc.readaxis = 'column'
    else:
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

    for science_file in science_list:
        iraf.noao.imred.ccdred.ccdproc(images = 'cleaned/cr.' + science_file,
                                       output = 'flattened/fl.' + science_file)
        #iraf.imcopy('cleaned/cr.' + science_file,'flattened/fl.' + science_file)

    iraf.unlearn(iraf.noao.twodspec.apextract.apall)

    if 'OSMOS' in instrument.upper():
        iraf.noao.twodspec.apextract.dispaxis = '1' # horizontal dispersion
    else:
        iraf.noao.twodspec.apextract.dispaxis = '2' # vertical dispersion
    iraf.noao.twodspec.apextract.apall.apertures = '1'
    iraf.noao.twodspec.apextract.apall.format = 'multispec'
    iraf.noao.twodspec.apextract.apall.references = ''
    iraf.noao.twodspec.apextract.apall.profiles = ''
    iraf.noao.twodspec.apextract.apall.interactive = apall_interactive
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
    iraf.noao.twodspec.apextract.apall.readnoise = readnoise
    iraf.noao.twodspec.apextract.apall.gain = gain
    iraf.noao.twodspec.apextract.apall.lsigma = 4
    iraf.noao.twodspec.apextract.apall.usigma = 4
    iraf.noao.twodspec.apextract.apall.nsubaps = '1'

    os.mkdir('extract')

    for science_file in science_list:
        iraf.noao.twodspec.apextract.apall(input = 'flattened/fl.' + science_file,
                                           output = 'extract/ex.' + science_file)

    # Extract a blue and red lamp using the first science object as our trace.

    if len(std_list)>0:
        trace_ref_image = std_list[0]
    else:
        trace_ref_image = science_list[0]

    iraf.noao.twodspec.apextract.apall.references = 'flattened/fl.' + trace_ref_image
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

    for lamp_file in lamp_list:
        iraf.noao.twodspec.apextract.apall(input = 'trimmed/tr.' + lamp_file,
                                           output = 'lamp.'+lamp_file)

    iraf.flprcache()
