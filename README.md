# Extinction-Maps-Correction

Here is the github respository for our paper: Validations and corrections of the SFD and Planck reddening maps based on LAMOST and Gaia data
Paper link:[]

## Supported dust maps
1. SFD 1998 [(Schlegel, Finkbeiner and Davis 1998)](http://doi.org/10.1086/305772)
2. Planck2014-Radiance and Tau [(Planck Collaboration et al. 2013)](http://doi.org/10.1051/0004-6361/201323195)
3. Planck2019-Tau [(Irfan et al. 2019)](http://doi.org/10.1051/0004-6361/201834394)

## User Guide

### Download the Respository

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

After dowlownd maps you want, please put them in the same directory with extinction_correct.py.

### Querying the original and corrected reddening values of Maps
Maps are queried using `astropy.coordinates.SkyCoord` objects.
