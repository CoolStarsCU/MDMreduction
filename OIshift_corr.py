import astropy.io.ascii as at
import astropy.io.fits as fits
import numpy as np
import matplotlib.pyplot as plt
import pdb

'''
This code was written with Python 3 in mind, but it works OK in Python 2.7, except for an input() vs. raw_input() issue.
'''

def fan(array, nfan=None, transpose=False):
    '''
    Turn a 1D array into a 2D array with nfan number of rows and len(array) number of columns.
    nfan - Number; number of rows.
    transpose - Boolean; whether to transpose 2D array.
    '''

    if nfan is None:
        nfan = len(array)

    fannedarray = np.tile(array, (nfan, 1))

    if transpose:
        outarray = fannedarray.T
    else:
        outarray = fannedarray

    return outarray


def fillarr(step, start, stop, fanned=None, transpose=False):
    '''
    Creates 1D array with values in the range [start-stop] in step steps (last element = stop when possible).

    fanned - Number; 1D array is turned into 2D array using the fan function.
    transpose - Boolean; whether the fan function transposes the output 2D array.
    '''

    if start > stop: return None
    if step == 0: return None

    filledarray = np.arange(start, stop + step, step)
    if filledarray[-1] > stop:
        filledarray = filledarray[:-1]

    if fanned is not None:
        outarray = fan(filledarray, fanned, transpose)
    else:
        outarray = filledarray

    return outarray


def ccpeak(arr1, arr2, radius=None, flagcf=False, flaglag=False):
    '''
    Locates precise location of the peak in the cross-correlation function between two vectors (arr1 and arr2). It is assumed that arr2 is the reference array.

    flagcf - Boolean; whether to return cf array.
    flaglag - Boolean; whether to return lag array.
    '''

    n = len(arr1)
    if radius is None: radius = 50

    # Normalize the arrays to get normalized cross correlation values
    arr1 = (arr1 - np.mean(arr1)) / (np.std(arr1) * len(arr1))
    arr2 = (arr2 - np.mean(arr2)) / np.std(arr2)

    lag = fillarr(1, -radius, radius)
    cf = np.zeros(len(lag))
    for i,tmplag in enumerate(lag):
        if i < 0:
            cf[i] = np.correlate(arr1[-tmplag:],arr2[:tmplag])[0]
        elif i == 0:
            cf[i] = np.correlate(arr1,arr2)[0]
        else:
            cf[i] = np.correlate(arr1[:-tmplag],arr2[tmplag:])[0]
    ind = np.argmax(cf)

    srad = 3
    sublag = lag[max(ind-srad, 0):min(ind+srad, 2*radius)+1]
    subcf = cf[max(ind-srad, 0):min(ind+srad, 2*radius)+1]

    a = np.polyfit(sublag, subcf, 2)
    a = a[::-1]

    maxlag = -a[1] / (2. * a[2])
    nlag = len(lag)
    if maxlag < lag[0]: maxlag = lag[0]
    if maxlag > lag[-1]: maxlag = lag[-1]

    outarray = [maxlag]
    if flagcf:
        outarray.append(cf)
    if flaglag:
        outarray.append(lag)

    return outarray



def main(ref_sky, input_file='to_shift.lis', output_file='OI_shifts.tbl'):
    '''
    Cross-correlates sky spectra of objects in input_file list against a chosen reference sky spectrum ref_sky, which it usually is a long exposure one.

    ref_sky - String; filename (without the .fits extension) of spectrum chosen as reference.
    input_file - String; filename of list of spectra to be shifted.
    output_file - String; filename of table with three columns: wavelength shift, 1-sigma error, and quality factor for each spectrum in input_file.
    '''

    from astropy.modeling import models, fitting

    plt.ion() # Prevents the pyplot figure window to block the terminal
    RED = '#b2182b'

    # Read list of spectra to be shifted
    infile = at.read(input_file, data_start=0)
    spectra = infile.columns[0].data
    obj_code = infile.columns[1].data
    name = infile.columns[2].data

    # Read reference sky spectrum
    skyfile = at.read('wavecal/sky.' + ref_sky, data_start=0)
    ref_wavelength = skyfile.columns[0].data
    ref_flux = skyfile.columns[1].data
    nskies = len(ref_flux)

    # Loop through all sky spectra
    outarray = []
    with open(output_file,"w") as f:
        for ispec,spec in enumerate(spectra):
            if obj_code[ispec] == 'lamp': continue

            # Read sky spectrum
            skyspec = at.read('wavecal/sky.' + spec, data_start=0)
            wavelength = skyspec.columns[0].data
            flux = skyspec.columns[1].data

            # Cross-correlate sky spectrum with reference sky spectrum
            pixel_shift, ccf = ccpeak(flux, ref_flux, 10, flagcf=True)
            max_ccf = np.max(ccf)
            median_ccf = np.median(ccf)
            quality_factor = (max_ccf - median_ccf) / median_ccf

            # Set up plot
            plt.close()
            fig = plt.figure(1, figsize=(9,7))
            plt.clf()
            ax = fig.add_subplot(111)

            # Plot sky spectrum
            ax.set_title(spec)
            irange = np.where((wavelength >= 5552.) & (wavelength <=5592.))[0]
            ax.plot(wavelength[irange], flux[irange], drawstyle='steps-mid',
                    color='k', linewidth=1.1)
            ax.axvline(5577.399, linestyle='--', color='k', alpha=0.8)

            if len(irange) > 2 and quality_factor > 0.075:
                # Fit Gaussian+Constant function to spectrum
                imax = np.argmax(flux[irange])
                Gaus_init = models.Gaussian1D(amplitude=np.max(flux[irange]),
                                              mean=wavelength[irange][imax],
                                              stddev=3.)
                Cons_init = models.Const1D(amplitude=np.median(flux[irange]))
                fitter = fitting.LevMarLSQFitter()
                f_init = Gaus_init + Cons_init # Actual function to fit
                func = fitter(f_init, wavelength[irange], flux[irange])

                # Read output of LSQ fitter
                OH_center = func.mean_0.value
                shift = 5577.34 - OH_center
                if fitter.fit_info['param_cov'] is None:
                    shift = 0
                    OH_center_err = 10.
                else:
                    diag_cov = np.diag(fitter.fit_info['param_cov'])
                    OH_center_err = np.sqrt(diag_cov[1])

                # Plot fit
                ax.plot(wavelength[irange]+shift, flux[irange], color=RED,
                        drawstyle='steps-mid', linestyle='--', linewidth=1.1)
                plt.pause(0.05)
            else:
                shift = 0.
                OH_center_err = 10.

            # Show user plot, shift and error values found from fit
            print(spec)
            print('Quality factor: ' + format(quality_factor, '.5f'))
            print('Shift and error: ', format(shift, '.5f'),
                  format(OH_center_err, '.5f'))
            plt.show()
            raw_input("Press [Enter] to continue.")
            # For Python 3, use input() instead of raw_input()
            plt.close()

            f.write("{:>14} {:>6} {:>18} {:.4f} {:.5f} {:.5f}\n".format(
                    spec, obj_code[ispec], name[ispec],
                    shift, OH_center_err, quality_factor))

        f.close()
        return None
