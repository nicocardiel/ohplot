# ohplot
The script ohplot.py overplots the expected OH emission lines and a
user-defined spectrum. The latter can be:

- any of the 34 SDSS spectral cross-correlation templates (see
  http://classic.sdss.org/dr7/algorithms/spectemplates/index.html). In this
  case the parameter "--sdss n" must be given, with n running from 0 to 33.

- any of the 12 templates from the Kinney-Calzetti spectral atlas of galaxies
  (see http://www.stsci.edu/hst/observatory/crds/cdbs_kc96.html). In this case
  the parameter "--kinn n" must be employed, with n running from 0 to 11.

- a particular spectrum provided by the user in ASCII format (two columns:
  wavelength, in Angstroms, and flux). In this case, the parameter "--ascii
  filename" must be used, being "filename" the name of the ASCII file.

Both the SDSS and the Kinney-Calzetti templates are stored in the data
subdirectory.

The template spectrum is normalized before plotting.

The displayed OH lines correspond to the Iraf line list and to the Oliva et al.
(2013; see https://ui.adsabs.harvard.edu/abs/2013A%26A...555A..78O/abstract)
list. Note that both lists are not identical.

Usage:

```bash
$ python ohplot.py --kinn 10 --redshift 1.6
$ python ohplot.py --sdss 27 --redshift 1.6 --emlines
$ python ohplot.py --ascii example.ascii --redshift 1.6
```

To see the list of additional parameters:

```bash
$ python ohplot.py --help
```


