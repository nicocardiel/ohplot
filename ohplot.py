# Copyright 2019-2026 Universidad Complutense de Madrid
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
from pathlib import Path
from scipy.stats import norm
import sys

WVMIN_DEFAULT = 8900
WVMAX_DEFAULT = 25000

NTEMP_SDSS = 33
kinney_list = [
    "kinn-00_elliptical_template.ascii",
    "kinn-01_bulge_template.ascii",
    "kinn-02_s0_template.ascii",
    "kinn-03_sa_template.ascii",
    "kinn-04_sb_template.ascii",
    "kinn-05_sc_template.ascii",
    "kinn-06_starb1_template.ascii",
    "kinn-07_starb2_template.ascii",
    "kinn-08_starb3_template.ascii",
    "kinn-09_starb4_template.ascii",
    "kinn-10_starb5_template.ascii",
    "kinn-11_starb6_template.ascii",
]
NTEMP_KINN = len(kinney_list)


def convolve_oh_lines(ohlines_wave, ohlines_flux, sigma, crpix1, crval1, cdelt1, naxis1):
    xwave = crval1 + (np.arange(naxis1) + 1 - crpix1) * cdelt1
    spectrum = np.zeros(naxis1)
    for wave, flux in zip(ohlines_wave, ohlines_flux):
        sp_tmp = gauss_box_model(x=xwave, amplitude=flux, mean=wave, stddev=sigma)
        spectrum += sp_tmp

    return xwave, spectrum


def gauss_box_model(x, amplitude=1.0, mean=0.0, stddev=1.0, hpix=0.5):
    """Integrate a Gaussian profile."""
    z = (x - mean) / stddev
    z2 = z + hpix / stddev
    z1 = z - hpix / stddev
    return amplitude * (norm.cdf(z2) - norm.cdf(z1))


def ohplot(
    ax=None,
    wvmin=None,
    wvmax=None,
    sdss=None,
    kinn=None,
    ascii_template=None,
    redshift=0.0,
    sigma=2.0,
    flux_afactor=0.0,
    flux_bfactor=1.0,
    emlines=False,
    title=None,
    noiraf=False,
    nodelta=False,
    emir_filters=False,
    echo=False,
):

    if echo:
        print("\033[1m\033[31mExecuting: " + " ".join(sys.argv) + "\033[0m\n")

    if ax is None:
        raise ValueError("Matplotlib axes 'ax' is None.")

    # ---
    # avoid incompatible options
    input1 = sdss is not None
    input2 = kinn is not None
    input3 = ascii_template is not None
    inputs = input1 + input2 + input3
    if inputs == 0:
        print("WARNING: no input template has been chosen")
    elif inputs > 1:
        print("ERROR: you can only choose an input template")

    # ---
    filename = None

    # ---
    if wvmin is None:
        wvmin = float(WVMIN_DEFAULT)

    if wvmax is None:
        wvmax = float(WVMAX_DEFAULT)

    # wavelength sampling
    crpix1 = 1.0
    crval1 = wvmin
    cdelt1 = 0.5
    naxis1 = int((wvmax - wvmin) / cdelt1) + 1

    # ---
    wave_template = None
    sp_template = None
    if sdss is not None:
        if 0 <= sdss <= NTEMP_SDSS - 1:
            # template spectrum
            filename = "data/spDR2-{:03d}.fit".format(sdss)
            with fits.open(filename, mode="readonly") as hdulist:
                template_header = hdulist[0].header
                template_data = hdulist[0].data
            # naxis1
            naxis1_template = template_header["NAXIS1"]
            # center wavelength (log10) or first pixel
            coeff0 = template_header["COEFF0"]
            # Log10 dispersion per pixel
            coeff1 = template_header["COEFF1"]
            # wavelength axis (logarithmic units)
            wave_template = 10 ** (coeff0 + np.arange(naxis1_template) * coeff1)
            # template spectrum
            sp_template = template_data[0, :]
        else:
            print("WARNING: template number out of range")
    elif kinn is not None:
        if 0 <= kinn <= NTEMP_KINN - 1:
            filename = "data/" + kinney_list[kinn]
            template_tabulated = np.genfromtxt(filename)
            wave_template = template_tabulated[:, 0]
            sp_template = template_tabulated[:, 1]
        else:
            print("WARNING: template number out of range")
    elif ascii_template is not None:
        filename = ascii_template
        template_tabulated = np.genfromtxt(ascii)
        wave_template = template_tabulated[:, 0]
        sp_template = template_tabulated[:, 1]

    # normalize flux to maximum value
    if sp_template is not None:
        sp_template /= sp_template.max()

    # ---
    # read sky OH lines from Oliva et al. (2013)
    ohlines_oliva = np.genfromtxt("data/Oliva_etal_2013.dat")

    # extract subset of lines within current wavelength range
    lok1 = ohlines_oliva[:, 1] >= wvmin
    lok2 = ohlines_oliva[:, 0] <= wvmax
    ohlines_oliva = ohlines_oliva[lok1 * lok2]

    # define wavelength and flux as separate arrays
    ohlines_oliva_wave = np.concatenate((ohlines_oliva[:, 1], ohlines_oliva[:, 0]))
    ohlines_oliva_flux = np.concatenate((ohlines_oliva[:, 2], ohlines_oliva[:, 2]))
    ohlines_oliva_flux /= ohlines_oliva_flux.max()

    # convolve location of OH lines to generate expected spectrum
    xwave_oliva, sp_oh_oliva = convolve_oh_lines(
        ohlines_oliva_wave, ohlines_oliva_flux, sigma, crpix1, crval1, cdelt1, naxis1
    )

    # normalize flux to maximum value
    sp_oh_oliva /= sp_oh_oliva.max()

    # ---
    # read sky OH lines from Iraf file
    ohlines_iraf = np.genfromtxt("data/ohlines_iraf_FULL.dat")

    # extract subset of lines within current wavelength range
    lok1 = ohlines_iraf[:, 0] >= wvmin
    lok2 = ohlines_iraf[:, 0] <= wvmax
    ohlines_iraf = ohlines_iraf[lok1 * lok2]

    # define wavelength and flux as separate arrays
    ohlines_iraf_wave = ohlines_iraf[:, 0]
    ohlines_iraf_flux = ohlines_iraf[:, 1]

    # convolve location of OH lines to generate expected spectrum
    xwave_iraf, sp_oh_iraf = convolve_oh_lines(
        ohlines_iraf_wave, ohlines_iraf_flux, sigma, crpix1, crval1, cdelt1, naxis1
    )

    # normalize flux to maximum value
    sp_oh_iraf /= sp_oh_iraf.max()

    # ---
    # read telluric transmission
    telluric_tabulated = np.genfromtxt("data/skycalc_transmission_R20000.txt")
    xtelluric = telluric_tabulated[:, 0] * 10  # convert from nm to Angstrom
    ytelluric = telluric_tabulated[:, 1]
    ytelluric /= ytelluric.max() * 0.7

    # ---
    # overplot telluric transmission
    ax.plot(xtelluric, ytelluric, color="gray", linestyle="-", label="telluric", alpha=0.5)

    # overplot OH lines from Oliva et al. (2003)
    if not nodelta:
        ax.stem(ohlines_oliva_wave, ohlines_oliva_flux, linefmt="C4-", markerfmt=" ", basefmt="C4-", label="OH (O2003)")

    # overplot convolved iraf lines
    if not noiraf:
        ax.plot(xwave_iraf, sp_oh_iraf, "C2-", label="OH (iraf)")

    # overplot convolved Oliva et al. (2003) lines
    ax.plot(xwave_oliva, sp_oh_oliva, "C1-", label="OH (O2003)")

    # EMIR filters
    if emir_filters:
        for filter_name, color in zip(["YJ", "HK", "K"], ["C0", "C2", "C3"]):
            emir_table = np.genfromtxt(f"data/filter_EMIR_{filter_name}spec.txt")
            wave_emir_filter = emir_table[:, 0] * 10000  # convert from microns to Angstroms
            transmission_emir_filter = emir_table[:, 1] / 100  # convert from percentage to fraction
            ax.plot(wave_emir_filter, transmission_emir_filter, color=color, label=f"EMIR {filter_name} filter")

    # plot limits
    ax.set_xlim([wvmin, wvmax])
    ymin, ymax = ax.get_ylim()
    dy = ymax - ymin
    ymin = 0
    ymax = ymax + 0.05 * dy
    ax.set_ylim(ymin, ymax)

    # typical emission lines (vacuum wavelengths in Angstroms)
    if emlines:
        emission_lines = [
            "[OII],                  3727.092",
            "[OII],                  3729.875",
            "${\\rm H}\\beta$,       4862.721",
            "[OIII],                 4960.295",
            "[OIII],                 5008.239",
            "[NII] ,                 6549.860",
            "${\\rm H}\\alpha$,      6564.614",
            "[NII],                  6585.270",
            "[SII],                  6718.290",
            "[SII],                  6732.680",
            "[FeII],                12570.21",
            "${\\rm Pa}\\beta$,     12821.69",
            "[FeII],                15330.0",
            "[FeII],                16440.0",
            "${\\rm H}_2$ S(5),     18345.0",
            "${\\ He}$I,            18635.0",
            "${\\rm Pa}\\alpha$,    18756.1",
            "[Si XI],               19320.0",
            "${\\rm Br}\\delta$,    19446.1",
            "${\\rm H}_2$ 1-0 S(3), 19576.0",
            "[Si VI],               19630.0",
            "${\\rm H}_2$ 1-0 S(2), 20338.0",
            "${\\rm H}_2$ 1-0 S(1), 21218.0",
            "${\\rm Br}\\gamma$,    21661.2",
            "${\\rm H}_2$ 1-0 S(0), 22235.0",
            "${\\rm H}_2$ 2-1 S(1), 22471.0",
        ]
        dy = ymax - ymin
        nplot = 0
        for item in emission_lines:
            namedum, wavedum = item.rsplit(",", maxsplit=1)
            wavedum = float(wavedum)
            if wvmin < wavedum * (1 + redshift) < wvmax:
                nplot += 1
                if nplot % 2 == 0:
                    delta_text = 0.04 * dy
                    color = "C0"
                else:
                    delta_text = 0.08 * dy
                    color = "black"
                xdum = wavedum * (1 + redshift)
                ax.plot([xdum, xdum], [0, ymax - 0.10 * dy], "--", linewidth=1, color=color)
                print(f"\nOverplotting emission line {namedum}")
                print(f"- Rest wavelength.... (Angstroms): {wavedum:.2f}")
                print(f"- Observed wavelength (Angstroms): {xdum:.2f}")
                ax.text(
                    wavedum * (1 + redshift),
                    ymax - delta_text,
                    rf"{namedum}",
                    fontsize=10,
                    horizontalalignment="center",
                    verticalalignment="bottom",
                    color=color
                )

    # overplot template spectrum
    if wave_template is not None:
        sp_scaled = flux_afactor + sp_template * flux_bfactor
        ax.plot(wave_template * (1 + redshift), sp_scaled, "w-", linewidth=3)
        ax.plot(wave_template * (1 + redshift), sp_scaled, "C7-", linewidth=1)

    # plot labels
    ax.set_xlabel("wavelength (Angstrom; in vacuum)")
    ax.set_ylabel("flux (arbitrary units)")

    # plot legend
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.15), ncol=4, fancybox=True, shadow=True)

    if filename is not None:
        ax.text(
            0.0,
            1.02,
            filename,
            fontsize=10,
            horizontalalignment="left",
            verticalalignment="bottom",
            transform=ax.transAxes,
        )

    ax.text(
        1.0,
        1.02,
        "z={:8.6f}".format(redshift),
        fontsize=10,
        horizontalalignment="right",
        verticalalignment="bottom",
        transform=ax.transAxes,
    )

    if title is not None:
        ax.set_title(title, fontsize=12)


def main(args=None):

    # parse command-line options
    parser = argparse.ArgumentParser(description="Overplot OH lines")

    # optional arguments
    parser.add_argument("--wvmin", help="Minimum wavelength", type=float)
    parser.add_argument("--wvmax", help="Maximum wavelength", type=float)
    parser.add_argument("--sdss", help="SDSS template number", type=int)
    parser.add_argument("--kinn", help="Kinney-Calzetti template number", type=int)
    parser.add_argument(
        "--ascii-template",
        help="ASCII file with template spectrum. Two columns: wavelength (Angstroms) and flux",
        type=str,
    )
    parser.add_argument("--redshift", help="Redshift for template spectrum", type=float, default=0.0)
    parser.add_argument(
        "--sigma", help="Broadening sigma (Angstroms) for OH lines (default: 2.0)", type=float, default=2.0
    )
    parser.add_argument("--flux_afactor", help="Additive factor for template spectrum", type=float, default=0.0)
    parser.add_argument("--flux_bfactor", help="Multiplicative factor for template spectrum", type=float, default=1.0)
    parser.add_argument("--title", help="Title for the plot", type=str)
    parser.add_argument("--emlines", help="overplot typical emission lines", action="store_true")
    parser.add_argument("--emir-filters", help="overplot EMIR filters", action="store_true")
    parser.add_argument("--nodelta", help="do not overplot delta functions for OH lines", action="store_true")
    parser.add_argument("--noiraf", help="do not overplot OH lines from Iraf", action="store_true")
    parser.add_argument(
        "--figsize", help="Figure size (width, height; default: 8.0 6.0)", type=float, nargs=2, default=[8.0, 6.0]
    )
    parser.add_argument("--echo", help="Display full command line", action="store_true")

    args = parser.parse_args()

    ascii_template = args.ascii_template
    if ascii_template is not None:
        if not Path(ascii_template).is_file():
            print(f"ERROR: ASCII file '{ascii_template}' does not exist")
            sys.exit(1)

    fig, ax = plt.subplots(figsize=(args.figsize[0], args.figsize[1]))
    ohplot(
        ax=ax,
        wvmin=args.wvmin,
        wvmax=args.wvmax,
        sdss=args.sdss,
        kinn=args.kinn,
        ascii_template=ascii_template,
        redshift=args.redshift,
        sigma=args.sigma,
        flux_afactor=args.flux_afactor,
        flux_bfactor=args.flux_bfactor,
        emlines=args.emlines,
        title=args.title,
        noiraf=args.noiraf,
        nodelta=args.nodelta,
        emir_filters=args.emir_filters,
        echo=args.echo,
    )
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":

    main()
