#!/usr/bin/env python
#!/user/covey/iraf/mypython

'''
This script contains functions that correct images for bias.
The output files have the an added "b" to their names.
'''

import os
import glob
import shutil
import numpy as np
import string as str
from sys import argv, exit
from pyraf import iraf
from pyraf.iraf import noao
from pyraf.iraf import imred, bias, colbias
import astropy.io.fits as pyfits
from list_utils import read_reduction_list
import pdb

def proc4k(argv='to_subtract.lis', instrument='OSMOSr4k', suffix='b', overwrite=False):
    '''
    This function is good for both MDM4K and MDMR4K images.

    This function is copied from proc4k.py by Paul Martini (OSU) and slightly edited.
    It performs the overscan subtraction and remove the relative gain differences for a single R4K image or a list of R4K images.
    Steps:
    1. determine if input is a file or list of files
    2. identify binning, size of overscan
    3. remove overscan and trim
    4. remove relative gain variations (Note: Gain for MDMR4K detector is still unknown as of 2017-11.)

    argv - String; it can be the name of a fits file, or the name of an ascii file that contains a list of fits filenames.
    instrument - String; it can be either OSMOS4k or OSMOSr4k.
    suffix - String; suffix to add to output filename(s).
    overwrite - Boolean; whether to overwrite the input file(s) with the output. Note that if overwrite=True and suffix is not None, the name of the overwritten file with be updated to reflect the suffix specified.

    Version 1.2: made it Python 3.0 friendly. Modified script to turn it into a function, with three input parameters. Parameter argv was modified to make it possible to input a single filename or a list of filenames in an ascii file. Read-out of fits image data is now explicit to be a numpy.int32 type.
    Version 1.3: made it flexible to be used with either OSMOS4k or OSMOSrk4 images. Added new instrument parameter.
    '''

    versNum = "1.3.0"
    versDate = "2017-12-11"

    # Define some variables
    Debug = False
    BiasSingle = 0
    BiasRow = 1
    BiasFit = 2
    #BiasType = BiasRow
    BiasType = BiasSingle
    Gain = False    # keep as False until gain values are known
    if instrument.upper() == 'OSMOSR4K':
        R4K = True
    elif instrument.upper() == 'OSMOS4K':
        R4K = False
    else:
        print('Instrument must be either OSMOS4k or OSMOSr4k.')
        return

    # Gain values for each amplifier [to be computed]
    r4k_gain_q1e = 1.0
    r4k_gain_q1o = 1.0
    r4k_gain_q2e = 1.0
    r4k_gain_q2o = 1.0
    r4k_gain_q3e = 1.0
    r4k_gain_q3o = 1.0
    r4k_gain_q4e = 1.0
    r4k_gain_q4o = 1.0
    mdm4k_gain_q1 = 1.0
    mdm4k_gain_q2 = 1.0
    mdm4k_gain_q3 = 1.0
    mdm4k_gain_q4 = 1.0

    # switch to more primitive (slower) code at MDM
    AT_MDM = False
    user = os.getlogin()
    if str.find(user, 'obs24m') >= 0 or str.find(user, 'obs13m') >= 0:
      AT_MDM = True

    # Determine type of input
    files = []
    if '.fits' in argv:
        files.append(argv)
    else:
        with open(argv, 'r') as infile:
            for line in infile:
                files.append(line[:-1]) # [:-1] is to remove the \n

    # Start procedure
    for fl in files:
      if os.path.isfile(fl):
        fitsfile = pyfits.open(fl)
        naxis1 = fitsfile[0].header['NAXIS1']
        naxis2 = fitsfile[0].header['NAXIS2']
        overscanx = fitsfile[0].header['OVERSCNX']
        overscany = fitsfile[0].header['OVERSCNY']	# should be 0
        ccdxbin = fitsfile[0].header['CCDXBIN']
        ccdybin = fitsfile[0].header['CCDYBIN']
        detector = fitsfile[0].header['DETECTOR']
        telescope = fitsfile[0].header['TELESCOP']
        overscanx /= ccdxbin
        overscany /= ccdybin
        # OSMOS or direct? [useful for knowing if MIS keywords have values]
        OSMOS = True
        if str.find(telescope, 'McGraw') >= 0:
          OSMOS = False	# direct image with the 1.3m
        print('Processing %s[%d:%d] OVERSCANX=%d OVERSCANY=%d from %s obtained at the %s' % (fl, naxis1, naxis2, overscanx, overscany, detector, telescope))
        if overscanx * ccdxbin < 32:
          print('Error: OVERSCNX=%d less than 32 in %s' % (overscanx, fl))
          exit(1)
        if overscany > 0:
          print('Error: code not tested with OVERSCNY > 0')
          exit(1)

        # This is now superfluous: -----------------
        #if str.find(detector, 'R4K') < 0:
        #  # if not R4K, assume MDM4K
        #  R4K  = False
        # ------------------------------------------

        #   IRAF units: 1:32, 33:556, 557:1080, 1081:1112
        # Python units: 0:31, 32:555, 556:1079, 1080:1111
        c1 = overscanx 		        # 32   first image column counting from *zero*
        c2 = int(0.5 * naxis1) - 1	# 555  last image column on first half
        c3 = c2+1			        # 556  first image column on second half
        c4 = naxis1 - overscanx - 1 # 1079 last image column
        r1 = overscany 		        # 0    first image row
        r2 = int(0.5 * naxis2) - 1	# 523  last image row on first half
        r3 = r2+1			        # 524  first image row on second half
        r4 = naxis2-overscany-1  	# 1047 last image row
        outnaxis1 = c4 - c1 + 1		# 1048 columns in output, trimmed image
        outnaxis2 = r4 - r1 + 1		# 1048 rows in output, trimmed image
        collen = int(0.5*outnaxis1)	# number of rows in an image quadrant
        rowlen = int(0.5*outnaxis2)	# number of rows in an image quadrant
        #
        # Assumed layout: (ds9 perspective)
        #
        #    q2	q4
        #
        #	 q1	q3
        #
        # each R4K quadrant has an even 'e' and an odd 'o' amplifier
        #

        if Debug:
          print(' Quadrants in IRAF pixels: ')
          print(' q1: [%d:%d,%d:%d]' % (c1+1, c2+1, r1+1, r2+1))
          print(' q2: [%d:%d,%d:%d]' % (c1+1, c2+1, r3+1, r4+1))
          print(' q3: [%d:%d,%d:%d]' % (c3+1, c4+1, r1+1, r2+1))
          print(' q4: [%d:%d,%d:%d]' % (c3+1, c4+1, r3+1, r4+1))
        # Calculate the bias level for each amplifier
        data = fitsfile[0].data.astype(np.int32) # A different Numpy version may read the data array as type uint16. Not good.
        # identify the columns to use to calculate the bias level
        # skip the first and last columns of the overscan
        # changed to 'list' for hiltner due to primitive python version
        starti = 4/ccdxbin
        if AT_MDM:
          if R4K:
            cols_over_q1e = list(np.arange(starti, overscanx-2, 2))
            cols_over_q1o = list(np.arange(starti+1, overscanx-2, 2))
            cols_over_q2e = cols_over_q1e
            cols_over_q2o = cols_over_q1o
            cols_over_q3e = list(np.arange(naxis1-overscanx+starti, naxis1-2, 2))
            cols_over_q3o = list(np.arange(naxis1-overscanx+starti+1, naxis1-2, 2))
            cols_over_q4e = cols_over_q3e
            cols_over_q4o = cols_over_q3o
            cols_q1e = list(np.arange(c1,c2,2))
            cols_q1o = list(np.arange(c1+1,c2+2,2))
            cols_q2e = cols_q1e
            cols_q2o = cols_q1o
            cols_q3e = list(np.arange(c3,c4,2))
            cols_q3o = list(np.arange(c3+1,c4+2,2))
            cols_q4e = cols_q3e
            cols_q4o = cols_q3o
          else:
            cols_over_q1 = list(np.arange(starti, overscanx-2, 1))
            cols_over_q2 = cols_over_q1
            cols_over_q3 = list(np.arange(naxis1-overscanx+starti, naxis1-2, 1))
            cols_over_q4 = cols_over_q3
            cols_q1 = list(np.arange(c1,c2+1,1))
            cols_q2 = cols_q1
            cols_q3 = list(np.arange(c3,c4+1,1))
            cols_q4 = cols_q3
        else:
          if R4K:
            # identify the even and odd columns in the overscan
            cols_over_q1e = np.arange(starti, overscanx-starti, 2)
            cols_over_q1o = np.arange(starti+1, overscanx-starti, 2)
            cols_over_q2e = cols_over_q1e
            cols_over_q2o = cols_over_q1o
            cols_over_q3e = np.arange(naxis1-overscanx+starti, naxis1-starti, 2)
            cols_over_q3o = np.arange(naxis1-overscanx+starti+1, naxis1-starti, 2)
            cols_over_q4e = cols_over_q3e
            cols_over_q4o = cols_over_q3o
            # identify the even and odd columns in each quadrant
            cols_q1e = np.arange(c1,c2,2)
            cols_q2e = cols_q1e
            cols_q1o = np.arange(c1+1,c2+2,2)
            cols_q2o = cols_q1o
            cols_q3e = np.arange(c3,c4,2)
            cols_q4e = cols_q3e
            cols_q3o = np.arange(c3+1,c4+2,2)
            cols_q4o = cols_q3o
          else:
            cols_over_q1 = np.arange(starti, overscanx-2, 1)
            cols_over_q2 = cols_over_q1
            cols_over_q3 = np.arange(naxis1-overscanx+starti, naxis1-2, 1)
            cols_over_q4 = cols_over_q3
            cols_q1 = np.arange(c1,c2+1,1)
            cols_q2 = cols_q1
            cols_q3 = np.arange(c3,c4+1,1)
            cols_q4 = cols_q3
        if Debug:
          print('Overscan columns: ')
          print('Q1/Q2 overscan even first and last columns:', (cols_over_q1e[0], cols_over_q1e[-1], len(cols_over_q1e)))
          print('Q1/Q2 overscan odd first and last columns:', (cols_over_q1o[0], cols_over_q1o[-1], len(cols_over_q1o)))
          print('Q3/Q4 overscan even first and last columns:', (cols_over_q3e[0], cols_over_q3e[-1], len(cols_over_q3e)))
          print('Q3/Q4 overscan odd first and last columns:', (cols_over_q3o[0], cols_over_q3o[-1], len(cols_over_q3o)))
        if Debug:
          print('Image columns: ')
          print('Q1/Q2 even first and last columns:', (cols_q1e[0], cols_q1e[-1], len(cols_q1e), r1, r2, len(cols_q1e)))
          print('Q1/Q2 odd first and last columns:', (cols_q1o[0], cols_q1o[-1], len(cols_q1o), r1+rowlen, r2+rowlen, len(cols_q1o)))
          print('Q3/Q4 even first and last columns:', (cols_q3e[0], cols_q3e[-1], len(cols_q3e), r1, r2, len(cols_q3e)))
          print('Q3/Q4 odd first and last columns:', (cols_q3o[0], cols_q3o[-1], len(cols_q3o), r1+rowlen, r2+rowlen, len(cols_q3o)))
        # create arrays with the median overscan vs. row for each amplifier
        if R4K:
          bias_q1e = np.zeros(rowlen, dtype=float)
          bias_q1o = np.zeros(rowlen, dtype=float)
          bias_q2e = np.zeros(rowlen, dtype=float)
          bias_q2o = np.zeros(rowlen, dtype=float)
          bias_q3e = np.zeros(rowlen, dtype=float)
          bias_q3o = np.zeros(rowlen, dtype=float)
          bias_q4e = np.zeros(rowlen, dtype=float)
          bias_q4o = np.zeros(rowlen, dtype=float)
        else:
          bias_q1 = np.zeros(rowlen, dtype=float)
          bias_q2 = np.zeros(rowlen, dtype=float)
          bias_q3 = np.zeros(rowlen, dtype=float)
          bias_q4 = np.zeros(rowlen, dtype=float)
        # calculate 1-D bias arrays for each amplifier
        for i in range(r1, r2+1, 1):
          if R4K:
            bias_q1e[i] = np.median(data[i,cols_over_q1e]) 	# data[rows, columns]
            bias_q1o[i] = np.median(data[i,cols_over_q1o])
            bias_q2e[i] = np.median(data[i+rowlen,cols_over_q2e])
            bias_q2o[i] = np.median(data[i+rowlen,cols_over_q2o])
            bias_q3e[i] = np.median(data[i,cols_over_q3e])
            bias_q3o[i] = np.median(data[i,cols_over_q3o])
            bias_q4e[i] = np.median(data[i+rowlen,cols_over_q4e])
            bias_q4o[i] = np.median(data[i+rowlen,cols_over_q4o])
          else:
            bias_q1[i] = np.median(data[i,cols_over_q1]) 	# data[rows, columns]
            bias_q2[i] = np.median(data[i+rowlen,cols_over_q2])
            bias_q3[i] = np.median(data[i,cols_over_q3])
            bias_q4[i] = np.median(data[i+rowlen,cols_over_q4])

        ##########################################################################
        # Subtract the bias from the output
        ##########################################################################
        if BiasType == BiasSingle:
          OverscanKeyValue = 'BiasSingle'
          # subtract a single bias value for each amplifier
          if R4K:
            bq1e = int(np.median(bias_q1e))
            bq1o = int(np.median(bias_q1o))
            bq2e = int(np.median(bias_q2e))
            bq2o = int(np.median(bias_q2o))
            bq3e = int(np.median(bias_q3e))
            bq3o = int(np.median(bias_q3o))
            bq4e = int(np.median(bias_q4e))
            bq4o = int(np.median(bias_q4o))
            if AT_MDM:
              for r in range(r1,r2+1):
                for c in cols_q1e:
                  data[r,c] -= bq1e
                for c in cols_q1o:
                  data[r,c] -= bq1o
                for c in cols_q2e:
                  data[r+rowlen,c] -= bq2e
                for c in cols_q2o:
                  data[r+rowlen,c] -= bq2o
                for c in cols_q3e:
                  data[r,c] -= bq3e
                for c in cols_q3o:
                  data[r,c] -= bq3o
                for c in cols_q4e:
                  data[r+rowlen,c] -= bq4e
                for c in cols_q4o:
                  data[r+rowlen,c] -= bq4o
            else:
              data[r1:r2+1,cols_q1e] -= bq1e
              data[r1:r2+1,cols_q1o] -= bq1o
              data[r3:r4+1,cols_q2e] -= bq2e
              data[r3:r4+1,cols_q2o] -= bq2o
              data[r1:r2+1,cols_q3e] -= bq3e
              data[r1:r2+1,cols_q3o] -= bq3o
              data[r3:r4+1,cols_q4e] -= bq4e
              data[r3:r4+1,cols_q4o] -= bq4o
          else:
            bq1 = int(np.median(bias_q1))
            bq2 = int(np.median(bias_q2))
            bq3 = int(np.median(bias_q3))
            bq4 = int(np.median(bias_q4))
            if AT_MDM:
                  for r in range(r1,r2+1):
                    for c in cols_q1:
                      data[r,c] -= bq1
                    for c in cols_q2:
                      data[r+rowlen,c] -= bq2
                    for c in cols_q3:
                      data[r,c] -= bq3
                    for c in cols_q4:
                      data[r+rowlen,c] -= bq4
            else:
              data[r1:r2+1,cols_q1] -= bq1
              data[r3:r4+1,cols_q2] -= bq2
              data[r1:r2+1,cols_q3] -= bq3
              data[r3:r4+1,cols_q4] -= bq4
        elif BiasType == BiasRow:
          # not implemented on Hiltner, for MDM4K, etc.
          print('Warning: This mode has not been fully tested')
          OverscanKeyValue = 'BiasRow'
          # subtract a bias value for each row of each amplifier
          #print r1, r2, len(bias_q1e)
          for i in range(r1, r2, 1):
            data[i,cols_q1e] -= bias_q1e[i]
            data[i,cols_q1o] -= bias_q1o[i]
            data[i+rowlen,cols_q2e] -= bias_q2e[i]
            data[i+rowlen,cols_q2o] -= bias_q2o[i]
            data[i,cols_q3e] -= bias_q3e[i]
            data[i,cols_q3o] -= bias_q3o[i]
            data[i+rowlen,cols_q4e] -= bias_q4e[i]
            data[i+rowlen,cols_q4o] -= bias_q4o[i]
        elif BiasType == BiasFit:
          OverscanKeyValue = 'BiasFit'
          print('Error: Have not implemented a fit to the bias yet. Please use BiasSingle')
          exit(1)
        else:
          print('Error: Bias subtraction type not parsed correctly')
          exit(1)

        ##########################################################################
        # Apply the gain correction  [not yet implemented]
        ##########################################################################

        if Gain:
          if R4K:
            if AT_MDM:
              for r in range(r1,r2+1):
                for c in cols_q1e:
                  data[r,c] -= r4k_gain_q1e
                for c in cols_q1o:
                  data[r,c] -= r4k_gain_q1o
                for c in cols_q2e:
                  data[r+rowlen,c] -= r4k_gain_q2e
                for c in cols_q2o:
                  data[r+rowlen,c] -= r4k_gain_q2o
                for c in cols_q2o:
                  data[r,c] -= r4k_gain_q3e
                for c in cols_q2o:
                  data[r,c] -= r4k_gain_q3o
                for c in cols_q2o:
                  data[r+rowlen,c] -= r4k_gain_q4e
                for c in cols_q2o:
                  data[r+rowlen,c] -= r4k_gain_q4o
            else:
              data[r1:r2,cols_q1e] /= r4k_gain_q1e
              data[r1:r2,cols_q1o] /= r4k_gain_q1o
              data[r3:r4,cols_q2e] /= r4k_gain_q2e
              data[r3:r4,cols_q2o] /= r4k_gain_q2o
              data[r1:r2,cols_q3e] /= r4k_gain_q3e
              data[r1:r2,cols_q3o] /= r4k_gain_q3o
              data[r3:r4,cols_q4e] /= r4k_gain_q4e
              data[r3:r4,cols_q4o] /= r4k_gain_q4o
          else:
            if AT_MDM:
              for r in range(r1,r2+1):
                for c in cols_q1:
                  data[r,c] /= mdm4k_gain_q1
                for c in cols_q2:
                  data[r+rowlen,c] /= mdm4k_gain_q2
                for c in cols_q2:
                  data[r,c] /= mdm4k_gain_q3
                for c in cols_q2:
                  data[r+rowlen,c] /= mdm4k_gain_q4
            else:
              data[r1:r2,cols_q1] /= mdm4k_gain_q1
              data[r3:r4,cols_q2] /= mdm4k_gain_q2
              data[r1:r2,cols_q3] /= mdm4k_gain_q3
              data[r3:r4,cols_q4] /= mdm4k_gain_q4


        ##########################################################################
        # Write the output file
        ##########################################################################

        fitsfile[0].data = data[r1:r4+1,c1:c4+1]
        OverscanKeyComment = 'Overscan by proc4k.py v%s (%s)' % (versNum, versDate)
        GainKeyValue = 'Relative'
        GainKeyComment = 'Gain removed by proc4k.py'
        #BiasKeyValue = '%s' % (versNum)
        #BiasKeyComment = 'Gain removed by proc4k.py'

        if OSMOS:
            # Prevent a Pyfits error when these are not assigned values
            try:
                fitsfile[0].header.remove('MISFILT', ignore_missing=True)
                fitsfile[0].header.remove('MISFLTID', ignore_missing=True)
                fitsfile[0].header.insert('GPROBEX', ('MISFLTID', -1, 'MIS Filter ID'))
                fitsfile[0].header.insert('GPROBEX', ('MISFILT', -1, 'MIS Filter Number'))
            except:
                if Debug:
                    print('Note: MISFILT and MISFLTID keywords not found')

            # Correct the wrong date format in OSMOS R4K images
            tmpdate = fitsfile[0].header['DATE-OBS'].split('T')
            if len(tmpdate) > 1:
                fitsfile[0].header['DATE-OBS'] = tmpdate[0]

        fitsfile[0].header['BIASPROC'] = (OverscanKeyValue, OverscanKeyComment)
        #fitsfile[0].header.update('BIASVER', BiasKeyValue, BiasKeyComment)
        if R4K:
          fitsfile[0].header['BIASQ1E'] = (bq1e, 'Bias subtracted from Q1E')
          fitsfile[0].header['BIASQ1O'] = (bq1o, 'Bias subtracted from Q1O')
          fitsfile[0].header['BIASQ2E'] = (bq2e, 'Bias subtracted from Q2E')
          fitsfile[0].header['BIASQ2O'] = (bq2o, 'Bias subtracted from Q2O')
          fitsfile[0].header['BIASQ3E'] = (bq3e, 'Bias subtracted from Q3E')
          fitsfile[0].header['BIASQ3O'] = (bq3o, 'Bias subtracted from Q3O')
          fitsfile[0].header['BIASQ4E'] = (bq4e, 'Bias subtracted from Q4E')
          fitsfile[0].header['BIASQ4O'] = (bq4o, 'Bias subtracted from Q4O')
        else:
          fitsfile[0].header['BIASQ1'] = (bq1, 'Bias subtracted from Q1')
          fitsfile[0].header['BIASQ2'] = (bq2, 'Bias subtracted from Q2')
          fitsfile[0].header['BIASQ3'] = (bq3, 'Bias subtracted from Q3')
          fitsfile[0].header['BIASQ4'] = (bq4, 'Bias subtracted from Q4')

        if Gain:
          if R4K:
            fitsfile[0].header['GAINPROC'] = (GainKeyValue, GainKeyComment)
            fitsfile[0].header['GainQ1'] = (r4k_gain_q1, 'Gain for Q1')
            fitsfile[0].header['GainQ2'] = (r4k_gain_q2, 'Gain for Q2')
            fitsfile[0].header['GainQ3'] = (r4k_gain_q3, 'Gain for Q3')
            fitsfile[0].header['GainQ4'] = (r4k_gain_q4, 'Gain for Q4')

        # Set up the output filename
        outfile = fl[:str.find(fl, '.fits')] + suffix + '.fits'

        if overwrite:
            print('  Warning: Overwriting and renaming pre-existing file %s to %s' % (fl, outfile))
            os.remove(fl)

        if os.path.isfile(outfile):
            print('  Warning: Overwriting pre-existing file %s' % outfile)
            os.remove(outfile)
        fitsfile.writeto(outfile)
        fitsfile.close()
