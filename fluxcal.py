#!/usr/bin/env python
#!/user/covey/iraf/mypython

###############################################################################
#The notes that were here actually go with wavecal

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

import astropy.io.ascii as at

def fluxcal(imagelist):

# read in a list of spectra as defined in point 1 above:

    image_list = at.read(imagelist,data_start=0,
                         names=["ccdno","type","target","image_region",
                                "bias_region","reference_lamp"])
    spectranum = len(image_list)
        
# now create important subsets: flats, lamps, biases, objects and std spectra

# first make lists for all the objects
    science_list = []
    science_names = []
    std_list = []
    std_names = []
    
# fill those lists with the names of objects

    for j,target_type in enumerate(image_list["type"]):
            if target_type == 'obj':
                science_list.append(image_list["ccdno"][j])
                science_names.append(image_list["target"][j])
            elif target_type == 'std':
                science_list.append(image_list["ccdno"][j])
                science_names.append(image_list["target"][j])
                std_list.append(image_list["ccdno"][j])
                std_names.append(image_list["target"][j])
            else:
                continue

# find the number of objects in each list

    numscience = len(science_list)
    numstds = len(std_list)
    
# now do flux calibration - start with setting airmass to middle of exposure

    iraf.astutil.setairmass.observatory = 'kpno'
    iraf.astutil.setairmass.intype = 'middle'
    iraf.astutil.setairmass.outtype = 'effective'
    iraf.astutil.setairmass.ra = 'RA'
    iraf.astutil.setairmass.dec = 'DEC'
    iraf.astutil.setairmass.equinox = 'EQUINOX'
    iraf.astutil.setairmass.st = 'ST'
    iraf.astutil.setairmass.ut = 'TIME-OBS'
    iraf.astutil.setairmass.date = 'DATE-OBS'
    iraf.astutil.setairmass.exposure = 'EXPTIME'
    iraf.astutil.setairmass.airmass = 'AIRMASS'
    iraf.astutil.setairmass.utmiddle = 'MIDUT'
    iraf.astutil.setairmass.scale = '750.0'
    iraf.astutil.setairmass.show = 'yes'
    iraf.astutil.setairmass.update = 'yes'
    iraf.astutil.setairmass.override = 'yes'

    for science_file in science_list:
        iraf.astutil.setairmass(images='wavecal/dc.' + science_file)

# run standard on our standard star observations

    iraf.noao.onedspec.standard.samestar = 'no' 
    # Frequently observed different stars, changing that value
    iraf.noao.onedspec.standard.beam_switch = 'no'
    iraf.noao.onedspec.standard.apertures = ''
    iraf.noao.onedspec.standard.bandwidth = '20'
    iraf.noao.onedspec.standard.bandsep = '30'
    iraf.noao.onedspec.standard.fnuzero = '3.68E-20'
    iraf.noao.onedspec.standard.extinction = 'onedstds$kpnoextinct.dat'
    iraf.noao.onedspec.standard.caldir = 'onedstds$spec50cal/'
    iraf.noao.onedspec.standard.observatory = 'KPNO'
    iraf.noao.onedspec.standard.interact = 'yes'
    iraf.noao.onedspec.standard.graphics = 'stdgraph'
    iraf.noao.onedspec.standard.cursor = ''

    for i,std_file in enumerate(std_list):
        this_std_name = std_names[i].split(".")[0]
        iraf.noao.onedspec.standard.star_name = this_std_name
        iraf.noao.onedspec.standard(input = 'wavecal/dc.' + std_file, 
                                    output = 'stdfile')
    
# now run sensfunc to get sensitivity functions out

    iraf.noao.onedspec.sensfunc.apertures = ''
    iraf.noao.onedspec.sensfunc.ignoreaps = 'yes'
    iraf.noao.onedspec.sensfunc.logfile = 'sensfunclog'
    iraf.noao.onedspec.sensfunc.extinction = 'onedstds$kpnoextinct.dat'
    iraf.noao.onedspec.sensfunc.observatory = 'KPNO'
    iraf.noao.onedspec.sensfunc.function = 'chebyshev'
    iraf.noao.onedspec.sensfunc.order = '3'
    iraf.noao.onedspec.sensfunc.interactive = 'yes'
    iraf.noao.onedspec.sensfunc.graphs = 'sr'
    iraf.noao.onedspec.sensfunc.marks = 'plus cross box'
    iraf.noao.onedspec.sensfunc.colors = '2 1 3 4'
    iraf.noao.onedspec.sensfunc.cursor = ''
    iraf.noao.onedspec.sensfunc.device = 'stdgraph'

    iraf.noao.onedspec.sensfunc(standards = 'stdfile', sensitivity = 'sens', newextinction = 'extinct.dat')

# then use those sensitivity functions to calibrate the data

    iraf.noao.onedspec.calibrate.extinct = 'yes'
    iraf.noao.onedspec.calibrate.flux = 'yes'
    iraf.noao.onedspec.calibrate.extinction = 'onedstds$kpnoextinct.dat'
    iraf.noao.onedspec.calibrate.observatory = 'KPNO'
    iraf.noao.onedspec.calibrate.ignoreaps = 'yes'
    iraf.noao.onedspec.calibrate.fnu = 'no'
#    iraf.noao.onedspec.calibrate.airmass = 'QAIRMASS'
#    iraf.noao.onedspec.calibrate.exptime = 'QEXPTIME'

    # find the max and min regions used by the sensitivity function
    sensfunc_log = at.read("sensfunclog",data_start=6)
    sensfunc_wave = sensfunc_log['col1']

    # Final flux calibration 
    # and trimming beyond the good flux calibration region
    os.mkdir('finals')

    for j,science_file in enumerate(science_list):
        iraf.noao.onedspec.calibrate(input = 'wavecal/dc.' + science_file,
                                     output = 'finals/' + science_names[j],
                                     sensitivity = 'sens')
        iraf.noao.onedspec.scopy(input='finals/' + science_names[j],
                                 output = 'finals/trim.' + science_names[j],
                                 w1=sensfunc_wave[0],w2=sensfunc_wave[-1],
                                 rebin='no')
    
    iraf.flprcache()
