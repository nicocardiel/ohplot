# Copyright 2019-2022 Universidad Complutense de Madrid
#
# Author: Nicolás Cardiel (cardiel@ucm.es)
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE
#
# This script overplots the expected OH emission lines and
# a user-defined spectrum
#

import argparse
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
import sys

WVMIN_DEFAULT = 8900
WVMAX_DEFAULT = 25000

NTEMP_SDSS = 33
kinney_list = ['kinn-00_elliptical_template.ascii',
               'kinn-01_bulge_template.ascii',
               'kinn-02_s0_template.ascii',
               'kinn-03_sa_template.ascii',
               'kinn-04_sb_template.ascii',
               'kinn-05_sc_template.ascii',
               'kinn-06_starb1_template.ascii',
               'kinn-07_starb2_template.ascii',
               'kinn-08_starb3_template.ascii',
               'kinn-09_starb4_template.ascii',
               'kinn-10_starb5_template.ascii',
               'kinn-11_starb6_template.ascii']
NTEMP_KINN = len(kinney_list)


def convolve_oh_lines(ohlines_wave, ohlines_flux, sigma,
                      crpix1, crval1, cdelt1, naxis1):
    xwave = crval1 + (np.arange(naxis1) + 1 - crpix1) * cdelt1
    spectrum = np.zeros(naxis1)
    for wave, flux in zip(ohlines_wave, ohlines_flux):
        sp_tmp = gauss_box_model(x=xwave, amplitude=flux, mean=wave,
                                 stddev=sigma)
        spectrum += sp_tmp

    return xwave, spectrum


def gauss_box_model(x, amplitude=1.0, mean=0.0, stddev=1.0, hpix=0.5):
    """Integrate a Gaussian profile."""
    z = (x - mean) / stddev
    z2 = z + hpix / stddev
    z1 = z - hpix / stddev
    return amplitude * (norm.cdf(z2) - norm.cdf(z1))


def main(args=None):

    # parse command-line options
    parser = argparse.ArgumentParser(
        description='Overplot OH lines'
    )

    # optional arguments
    parser.add_argument("--wvmin",
                        help="Minimum wavelength",
                        type=float)
    parser.add_argument("--wvmax",
                        help="Maximum wavelength",
                        type=float)
    parser.add_argument("--sdss",
                        help="SDSS template number",
                        type=int)
    parser.add_argument("--kinn",
                        help="Kinney-Calzetti template number",
                        type=int)
    parser.add_argument("--ascii",
                        help="ASCII file with template spectrum",
                        type=argparse.FileType('rt'))
    parser.add_argument("--redshift",
                        help="Redshift for template spectrum",
                        type=float, default=0.0)
    parser.add_argument("--sigma",
                        help="Broadening sigma (Angstroms) for OH lines",
                        type=float, default=2.0)
    parser.add_argument("--flux_afactor",
                        help="Additive factor for template spectrum",
                        type=float, default=0.0)
    parser.add_argument("--flux_bfactor",
                        help="Multiplicative factor for template spectrum",
                        type=float, default=1.0)
    parser.add_argument("--emlines",
                        help="overplot typical emission lines",
                        action="store_true")
    parser.add_argument("--noiraf",
                        help="do not overplot OH lines from Iraf",
                        action="store_true")
    parser.add_argument("--echo",
                        help="Display full command line",
                        action="store_true")

    args = parser.parse_args()

    if args.echo:
        print('\033[1m\033[31mExecuting: ' + ' '.join(sys.argv) + '\033[0m\n')

    # ---
    # avoid incompatible options
    input1 = (args.sdss is not None)
    input2 = (args.kinn is not None)
    input3 = (args.ascii is not None)
    inputs = input1 + input2 + input3
    if inputs == 0:
        print('WARNING: no input template has been chosen')
    elif inputs > 1:
        print('ERROR: you can only choose an input template')

    # ---
    filename = None

    # ---
    if args.wvmin is None:
        wvmin = float(WVMIN_DEFAULT)
    else:
        wvmin = float(args.wvmin)

    if args.wvmax is None:
        wvmax = float(WVMAX_DEFAULT)
    else:
        wvmax = float(args.wvmax)

    # wavelength sampling
    crpix1 = 1.0
    crval1 = wvmin
    cdelt1 = 0.5
    naxis1 = int((wvmax-wvmin)/cdelt1) + 1

    # ---
    wave_template = None
    sp_template = None
    if args.sdss is not None:
        if 0 <= args.sdss <= NTEMP_SDSS:
            # template spectrum
            filename = 'data/spDR2-{:03d}.fit'.format(args.sdss)
            with fits.open(filename, mode='readonly') as hdulist:
                template_header = hdulist[0].header
                template_data = hdulist[0].data
            # naxis1
            naxis1_template = template_header['NAXIS1']
            # center wavelength (log10) or first pixel
            coeff0 = template_header['COEFF0']
            # Log10 dispersion per pixel
            coeff1 = template_header['COEFF1']
            # wavelength axis (logarithmic units)
            wave_template = 10 ** (coeff0 +
                                   np.arange(naxis1_template) * coeff1)
            # template spectrum
            sp_template = template_data[0, :]
        else:
            print('WARNING: template number out of range')
    elif args.kinn is not None:
        if 0 <= args.kinn <= NTEMP_KINN:
            filename = 'data/' + kinney_list[args.kinn]
            template_tabulated = np.genfromtxt(filename)
            wave_template = template_tabulated[:, 0]
            sp_template = template_tabulated[:, 1]
        else:
            print('WARNING: template number out of range')
    elif args.ascii is not None:
        filename = args.ascii.name
        template_tabulated = np.genfromtxt(args.ascii)
        wave_template = template_tabulated[:, 0]
        sp_template = template_tabulated[:, 1]

    # normalize flux to maximum value
    if sp_template is not None:
        sp_template /= sp_template.max()

    # ---
    # read sky OH lines from Oliva et al. (2013)
    ohlines_oliva = np.genfromtxt('data/Oliva_etal_2013.dat')

    # extract subset of lines within current wavelength range
    lok1 = ohlines_oliva[:, 1] >= wvmin
    lok2 = ohlines_oliva[:, 0] <= wvmax
    ohlines_oliva = ohlines_oliva[lok1 * lok2]

    # define wavelength and flux as separate arrays
    ohlines_oliva_wave = np.concatenate(
        (ohlines_oliva[:, 1], ohlines_oliva[:, 0]))
    ohlines_oliva_flux = np.concatenate(
        (ohlines_oliva[:, 2], ohlines_oliva[:, 2]))
    ohlines_oliva_flux /= ohlines_oliva_flux.max()

    # convolve location of OH lines to generate expected spectrum
    xwave_oliva, sp_oh_oliva = convolve_oh_lines(
        ohlines_oliva_wave, ohlines_oliva_flux, args.sigma,
        crpix1, crval1, cdelt1, naxis1
    )

    # normalize flux to maximum value
    sp_oh_oliva /= sp_oh_oliva.max()

    # ---
    # read sky OH lines from Iraf file
    ohlines_iraf = np.genfromtxt('data/ohlines_iraf_FULL.dat')

    # extract subset of lines within current wavelength range
    lok1 = ohlines_iraf[:, 0] >= wvmin
    lok2 = ohlines_iraf[:, 0] <= wvmax
    ohlines_iraf = ohlines_iraf[lok1 * lok2]

    # define wavelength and flux as separate arrays
    ohlines_iraf_wave = ohlines_iraf[:, 0]
    ohlines_iraf_flux = ohlines_iraf[:, 1]

    # convolve location of OH lines to generate expected spectrum
    xwave_iraf, sp_oh_iraf = convolve_oh_lines(
        ohlines_iraf_wave, ohlines_iraf_flux, args.sigma,
        crpix1, crval1, cdelt1, naxis1
    )

    # normalize flux to maximum value
    sp_oh_iraf /= sp_oh_iraf.max()

    # ---
    # read telluric transmission
    telluric_tabulated = np.genfromtxt('data/skycalc_transmission_R20000.txt')
    xtelluric = telluric_tabulated[:, 0] * 10  # convert from nm to Angstrom
    ytelluric = telluric_tabulated[:, 1]
    ytelluric /= ytelluric.max() * 0.7

    # ---
    fig, ax = plt.subplots()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width, box.height*0.95])

    # overplot telluric transmission
    ax.plot(xtelluric, ytelluric, 'C3-', label='telluric')

    # overplot OH lines from Oliva et al. (2003)
    ax.stem(ohlines_oliva_wave, ohlines_oliva_flux, linefmt='C4-',
            markerfmt=' ', basefmt='C4-', label='OH (O2003)')

    # overplot convolved iraf lines
    if not args.noiraf:
        ax.plot(xwave_iraf, sp_oh_iraf, 'C2-', label='OH (iraf)')

    # overplot convolved Oliva et al. (2003) lines
    ax.plot(xwave_oliva, sp_oh_oliva, 'C1-', label='OH (O2003)')

    # typical emission lines
    if args.emlines:
        emission_lines = [3727.092, 3729.875, 4862.721, 4960.295, 5008.239,
                          6549.860, 6564.614, 6585.270, 6718.290, 6732.680]
        for wavedum in emission_lines:
            ax.axvline(wavedum * (1 + args.redshift), color='k',
                       linestyle='--')

    # overplot template spectrum
    if wave_template is not None:
        sp_scaled = args.flux_afactor + sp_template * args.flux_bfactor
        ax.plot(wave_template * (1 + args.redshift), sp_scaled, 'C7-')

    # plot limits
    ax.set_xlim([wvmin, wvmax])

    # plot labels
    ax.set_xlabel('wavelength (Angstrom; in vacuum)')
    ax.set_ylabel('flux (arbitrary units)')

    # plot legend
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.20),
              ncol=4, fancybox=True, shadow=True)

    if filename is not None:
        ax.text(0.0, 1.02, filename, fontsize=10,
                horizontalalignment='left',
                verticalalignment='bottom',
                transform=ax.transAxes)
        ax.text(1.0, 1.02, 'z={:8.6f}'.format(args.redshift), fontsize=10,
                horizontalalignment='right',
                verticalalignment='bottom',
                transform=ax.transAxes)

    plt.show()


if __name__ == "__main__":

    main()
