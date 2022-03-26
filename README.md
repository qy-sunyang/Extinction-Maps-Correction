# Extinction-Maps-Correction

Here is the github respository for our paper: Validations and corrections of the SFD and Planck reddening maps based on LAMOST and Gaia data
Paper link:[]

## Supported dust maps
1. SFD 1998 [(Schlegel, Finkbeiner and Davis 1998)](http://doi.org/10.1086/305772)
2. Planck2014-Radiance and Tau [(Planck Collaboration et al. 2013)](http://doi.org/10.1051/0004-6361/201323195)
3. Planck2019-Tau [(Irfan et al. 2019)](http://doi.org/10.1051/0004-6361/201834394)

## User Guide

### Pre-requested Packages
1. [healpy (1.15.2)](https://healpy.readthedocs.io/en/latest/)
2. [astropy (4.3.1)](https://www.astropy.org/)

### Download the Whole Respository

### Download Original Maps
Here are links for downloading original maps:
1. SFD 1998 [https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/EWCNL5](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/EWCNL5)

Maps you will need are:
* SFD_dust_4096_ngp.fits 
* SFD_dust_4096_sgp.fits
* SFD_temp_ngp.fits
* SFD_temp_sgp.fits
2. Planck2014 [https://irsa.ipac.caltech.edu/data/Planck/release_1/all-sky-maps/previews/HFI_CompMap_ThermalDustModel_2048_R1.20/index.html](https://irsa.ipac.caltech.edu/data/Planck/release_1/all-sky-maps/previews/HFI_CompMap_ThermalDustModel_2048_R1.20/index.html)
3. Planck2019 [https://vizier.cds.unistra.fr/viz-bin/VizieR-4](https://vizier.cds.unistra.fr/viz-bin/VizieR-4)

Maps you will need are:
* temp.fits 
* beta.fits
* tau.fits

After dowlownd maps you want, please put them in the *data* directory.

### Querying the original and corrected reddening values of Maps
Maps are queried using `astropy.coordinates.SkyCoord` objects. Currently we only support 'icrs' and 'galactic' systems. For example,
```
>>> from astropy.coordinates import SkyCoord
>>> c = SkyCoord([180, 45., 180., 80.], [40,15, 0., -15.], frame="galactic", unit="deg")
>>> c = SkyCoord([180, 45., 180., 80.], [40,15, 0., -15.], frame="icrs", unit="deg")
```
Then call the function you need to obtain the original reddening values and corrected values. For example,
```
>>> from extinction_correct import sfd_reddening_correction
>>> from astropy.coordinates import SkyCoord
>>> c = SkyCoord([180, 45., 180., 80.], [40,15, 0., -15.], frame="galactic", unit="deg")
>>> out,out_correct = sfd_reddening_correction(c)
>>> print(out,out_correct)
[0.02857005 0.27206674 1.4048629  0.18249302] [0.02582535        nan        nan 0.15870613]
```
You may get warnings and nan for corrected results. Unfortunately, it means for specific coordinates we cannot provide corrected values.

If the warning infomation is *'Coordinates [i] are out of application range (E(B-V)>=0.3)'*, it means that,
1. Coordinates [i] are not within the footprint of LAMOST so that we cannot provide the correction factor k for those regions;
2. Since the original excess E(B-V) of Coordinates [i] are larger than 0.3, which is not at the application range of fitted k relations, we cannot provide the fitted correction factor k_fit as well.

If the warning infomation is *'Coordinates [1] cannot get a reasonable k_fit value.'*, it means that,
1. Coordinates [i] are not within the footprint of LAMOST so that we cannot provide the correction factor k for those regions;
2. They are at the application range of fitted k relations, but k_fit is out of range [0, 2] that is not reliable.
