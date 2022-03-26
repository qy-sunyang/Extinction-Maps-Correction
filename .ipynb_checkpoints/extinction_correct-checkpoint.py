import pandas as pd
import numpy as np
import healpy as hp
import os
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy import wcs
from scipy.ndimage import map_coordinates
from scipy import interpolate
import logging

def kfit_func(X,para):
    '''
    This function is foth-order polynomial function for k_fit.
    
    Args:
        1. X('list'): The variables in the k_fit functions.
        
        For SFD and Planck2014-R, X[0]('numpy.array') is the dust temperature 
        and X[1]('numpy.array') is the reddening value of the map for correction;
        
        For Planck2014-Tau and Planck2019-Tau, X[0]('numpy.ndarray') is the dust temperature 
        and X[1]('numpy.ndarray') is the spectral index(beta) of the map for correction.
        
        2. para('numpy.ndarray'): All parameters of the k_fit function. 
        
    Returns:
        1. out('numpy.float64'): Original reddening value E(B-V) from the SFD map;
            
        2. out_correct('numpy.float64'): The corrected reddening value from maps. 
        
    Warnings:
        1. 'Coordinates [i] are out of application range (E(B-V)>=0.3)': This 
        means this pixel is neither within the footprint of LAMOST where can 
        use k map correct the reddening directly, or within the proper range 
        of using the k_fit relationship.
        
        2. 'Coordinates [1] cannot get a reasonable k_fit value.': This means 
        the final k_fit of this pxiel is out of [0, 2], which we think it is 
        not a reasonable value for correction.
    '''
    
    x,y = X
    a, b, c, d, e, f, g, h, i, j, k, l, m, n, o = para
    return a*x**4 + b*y**4 + c*y*x**3 + d*x*y**3 + e*(x**2)*(y**2) + \
            f*y**3 + g*x**3 + h*y*x**2 + i*x*y**2 + j*x**2 +\
            k*y**2 + l*x*y + m*x + n*y + o

def sfd_reddening_correction(coordinate):
    '''
    This function is for correcting SFD extinction map.
    
    Args:
        1. coordinate('astropy.coordinates.SkyCoord'): The coordinates of target to correct extinction. 
            only two options can be accepted : 'icrs' or 'galactic'.
        
    Returns:
        1. out('numpy.float64'): Original reddening value E(B-V) from the SFD map;
            
        2. out_correct('numpy.float64'): The corrected reddening value from maps. 
        
    Warnings:
        1. 'Coordinates [i] are out of application range (E(B-V)>=0.3)': This means this pixel is neither 
        within the footprint of LAMOST where can use k map correct the reddening directly, or within the 
        proper range of using the k_fit relationship.
        
        2. 'Coordinates [1] cannot get a reasonable k_fit value.': This means the final k_fit of this pxiel
        is out of [0, 2], which we think it is not a reasonable value for correction.
    '''
    
    
    # convert coordinates to galactic frame     
    if coordinate.frame.name == 'icrs':
        coordinate = coordinate.galactic
        l = coordinate.galactic.l.deg
        b = coordinate.galactic.b.deg
    else:
        l = coordinate.l.deg
        b = coordinate.b.deg
    
    # prepare the origianl reddening value from the SFD map 
    pole = ['ngp','sgp']
    hduN=fits.open('data/SFD_dust_4096_ngp.fits')
    hduS=fits.open('data/SFD_dust_4096_sgp.fits')
    dataN = [hduN[0].data, wcs.WCS(hduN[0].header)]
    dataS = [hduS[0].data, wcs.WCS(hduS[0].header)]
    data_pole = [dataN,dataS]        
    out = np.full(len(l), np.nan, dtype='f4')
    for i in range(len(pole)):
        m = (b >= 0) if pole[i] == 'ngp' else (b < 0)
        if np.any(m):
            data, w = data_pole[i] 
            x, y = w.wcs_world2pix(l[m], b[m], 0)
            out[m] = map_coordinates(data, [y, x], order=1, mode='nearest')
            
    # Planck et al. 2013 indicates that SFD have an offset of −0.003 mag for the SFD map.           
    out_correct = out + 0.003 

    #load the k map of SFD     
    use_k_map = hp.fitsfunc.read_map\
    ('data/k_maps/k_SFD.fits', nest=False, dtype=None) 
    
    # convert coordinate to HEALPix indice  for k map    
    nside = 64 
    npix = hp.nside2npix(nside)
    theta = np.radians(90. - b)
    phi = np.radians(l)
    indices = hp.ang2pix(nside, theta, phi)
    
    ### read the value from k map ###
    k = use_k_map[indices]

    # check whether k is in the footprint. If not, k will be obtained from fitted k relation
    valid = (k > -1e+30)
    if np.any(valid):
        out_correct[valid] = out_correct[valid]*k[valid] 
        
    un_valid = (k < -1e+30)    
    if np.any(un_valid):
        unvalid_index = np.where(un_valid==True)[0]
        
        #check whether the coordinate is in the applicable range of EBV
        judge_outrange = unvalid_index[(out_correct[unvalid_index]>=0.3)]
        out_correct[judge_outrange]=np.nan
        logging.warning('Coordinates '+str(judge_outrange)+' are out of application range (E(B-V)>=0.3)')
        
        unvalid_index = unvalid_index[(out_correct[unvalid_index]<0.3)]
        l_forfit = l[unvalid_index]
        b_forfit = b[unvalid_index]
        
        # prepare SFD dust temperature map
        hduN=fits.open('data/SFD_temp_ngp.fits')
        hduS=fits.open('data/SFD_temp_sgp.fits')
        dataN = [hduN[0].data, wcs.WCS(hduN[0].header)]
        dataS = [hduS[0].data, wcs.WCS(hduS[0].header)]
        data_pole = [dataN,dataS]                
        temp = np.full(len(l_forfit), np.nan, dtype='f4')
        for i in range(len(pole)):
            m = (b_forfit >= 0) if pole[i] == 'ngp' else (b_forfit < 0)
            if np.any(m):
                data, w = data_pole[i] 
                x, y = w.wcs_world2pix(l_forfit[m], b_forfit[m], 0)
                temp[m] = map_coordinates(data, [y, x], order=1, mode='nearest')
                
        ### using SFD fit function to get correction factor k ###
        
        # read the values of fitting parameter:
        func = pd.read_csv('data/kfit_relation/all_func.csv')
        SFD_func = func['SFD'].values
        k_fit = kfit_func([temp,out_correct[unvalid_index]],SFD_func)#.values)
        out_correct[unvalid_index] = out_correct[unvalid_index] * k_fit
        
        # check whether k_fit is reasonable value 
        # if kfit is negative or too large. it cannot be used.        
        judge_kfit = unvalid_index[((k_fit<0)|(k_fit>2))] 
        if len(judge_kfit)!=0:
            out_correct[judge_kfit]=np.nan
            logging.warning('Coordinates '+str(judge_kfit)+' cannot get a reasonable k_fit value.')
    
    return out,out_correct


def Planck2014R_reddening_correction(coordinate):
    '''
    This function is for correcting Planck2014-R extinction map. (Planck Collaboration et al. 2013)
    
    Args:
        1. coordinate('astropy.coordinates.SkyCoord'): The coordinates of target to correct extinction. 
        Only two options can be accepted : 'icrs' or 'galactic'.
        
    Returns:
        1. out('numpy.float64'): Original reddening value E(B-V) from the Planck2014-R map;
            
        2. out_correct('numpy.float64'): The corrected reddening value from maps. 
        
    Warnings:
        1. 'Coordinates [i] are out of application range (E(B-V)>=0.3)': This means this pixel is neither 
        within the footprint of LAMOST where can use k map correct the reddening directly, or within the 
        proper range of using the k_fit relationship.
        
        2. 'Coordinates [i] cannot get a reasonable k_fit value.': This means the final k_fit of this pxiel
        is out of [0, 2], which we think it is not a reasonable value for correction.
    '''
    
    # convert coordinates to galactic frame 
    if coordinate.frame.name == 'icrs':
        coordinate = coordinate.galactic
        l = coordinate.galactic.l.deg
        b = coordinate.galactic.b.deg
    else:
        l = coordinate.l.deg
        b = coordinate.b.deg
    
    ### prepare the origianl reddening value from the Planck2014 map ###
    
    # index for field: TAU=0, RADIANCE=3, TEMP=4, BETA=6
    planck_map = hp.fitsfunc.read_map('data/HFI_CompMap_ThermalDustModel_2048_R1.20.fits',\
                                  field=3,nest=False,hdu=1)
    out = np.full(len(l), np.nan, dtype='f4')
    out_correct = np.full(len(l), np.nan, dtype='f4')  
    
    # convert Skycoord to indice in Planck2014 HEALPix map
    nside_planck = 2048
    theta = np.radians(90. - b)
    phi = np.radians(l)
    indices = hp.ang2pix(nside_planck, theta, phi)

    # 5.4*10**5 is the factor of the transformation between EBV and Radiance
    out = planck_map[indices]*5.4*10**5

        
    #load the k map of Planck2014-R 
    use_k_map = hp.fitsfunc.read_map\
    ('data/k_maps/k_Planck2014-R.fits', nest=False, dtype=None) 
    
    # convert Skycoord to indice in k HEALPix map     
    nside = 64 
    npix = hp.nside2npix(nside)
    theta = np.radians(90. - b)
    phi = np.radians(l)
    indices = hp.ang2pix(nside, theta, phi)
    
    ### read the value from k map ###
    k = use_k_map[indices]

    # check whether k is in the footprint. If not, k will be obtained from fitted k relation
    valid = (k > -1e+30)
    if np.any(valid):
        out_correct[valid] = out[valid]*k[valid]
        
    un_valid = (k < -1e+30)
    if np.any(un_valid):
        unvalid_index = np.where(un_valid==True)[0]
        
        # check whether the coordinate is in the applicable range of EBV
        judge_outrange = unvalid_index[(out[unvalid_index]>=0.3)]
        out_correct[judge_outrange]=np.nan
        logging.warning('Coordinates '+str(judge_outrange)+' are out of application range (E(B-V)>=0.3)')
        
        unvalid_index = unvalid_index[(out[unvalid_index]<0.3)]
        l_forfit = l[unvalid_index]
        b_forfit = b[unvalid_index]
        
        # Read Planck2014-R dust temperature map
        # TAU:0,RADIANCE:3,TEMP:4,BETA:6
        planck_temp = hp.fitsfunc.read_map('data/HFI_CompMap_ThermalDustModel_2048_R1.20.fits',\
                                  field=4,nest=False,hdu=1)
        temp = planck_temp[indices[unvalid_index]]

        ### using Planck2014-R fit function to get correction factor k ###
        
        # read the values of fitting parameter:
        func = pd.read_csv('data/kfit_relation/all_func.csv')
        Plank2014R_func = func['Planck2014-R'].values

        k_fit = kfit_func([temp,out[unvalid_index]],Plank2014R_func)
        out_correct[unvalid_index] = out[unvalid_index] * k_fit
        
        #check whether k_fit is reasonable value
        # if kfit is negative or too large. it cannot be used.
        judge_kfit = unvalid_index[((k_fit<0)|(k_fit>2))]        
        if len(judge_kfit)!=0:
            out_correct[judge_kfit]=np.nan
            logging.warning('Coordinates '+str(judge_kfit)+' cannot get a reasonable k_fit value.')

    return out,out_correct


def Planck2014Tau_reddening_correction(coordinate):
    '''
    This function is for correcting Planck2014-Tau extinction map. (Planck Collaboration et al. 2013)
    
    Args:
        1. coordinate('astropy.coordinates.SkyCoord'): The coordinates of target to correct extinction. 
        Only two options can be accepted : 'icrs' or 'galactic'.
        
    Returns:
        1. out('numpy.float64'): Original reddening value E(B-V) from the Planck2014-Tau map;
            
        2. out_correct('numpy.float64'): The corrected reddening value from maps. 
        
    Warnings:
        1. 'Coordinates [i] are out of application range (E(B-V)>=0.3)': This means this pixel is neither 
        within the footprint of LAMOST where can use k map correct the reddening directly, or within the 
        proper range of using the k_fit relationship.
        
        2. 'Coordinates [1] cannot get a reasonable k_fit value.': This means the final k_fit of this pxiel
        is out of [0, 2], which we think it is not a reasonable value for correction.
    '''
    
    # convert coordinates to galactic frame    
    if coordinate.frame.name == 'icrs':
        coordinate = coordinate.galactic
        l = coordinate.galactic.l.deg
        b = coordinate.galactic.b.deg
    else:
        l = coordinate.l.deg
        b = coordinate.b.deg
    
    # prepare the origianl reddening value from the map 
    # index for field: TAU=0, RADIANCE=3, TEMP=4, BETA=6
    planck_map = hp.fitsfunc.read_map('data/HFI_CompMap_ThermalDustModel_2048_R1.20.fits',\
                                  field=0,nest=False,hdu=1)
    out = np.full(len(l), np.nan, dtype='f4')
    out_correct = np.full(len(l), np.nan, dtype='f4')
    
    # convert Skycoord to indice in Planck2014 HEALPix map
    nside_planck = 2048
    theta = np.radians(90. - b)
    phi = np.radians(l)
    indices = hp.ang2pix(nside_planck, theta, phi)

    # 1.49*10**4 is the factor of the transformation between EBV and Tau
    out = planck_map[indices]*1.49*10**4

        
    #load the k map of Planck2014-Tau 
    use_k_map = hp.fitsfunc.read_map\
    ('data/k_maps/k_Planck2014-Tau.fits', nest=False, dtype=None) 
    
    # convert Skycoord to indice in k HEALPix map     
    nside = 64
    npix = hp.nside2npix(nside)
    theta = np.radians(90. - b)
    phi = np.radians(l)
    indices = hp.ang2pix(nside, theta, phi)
    
    ### read the value from k map ###
    k = use_k_map[indices]

    # check whether k is in the footprint. If not, k will be obtained from fitted k relation
    valid = (k > -1e+30)
    if np.any(valid):
        out_correct[valid] = out[valid]*k[valid]
        
    un_valid = (k < -1e+30)
    if np.any(un_valid):
        unvalid_index = np.where(un_valid==True)[0]
        
        # check whether the coordinate is in the applicable range of EBV
        judge_outrange = unvalid_index[(out[unvalid_index]>=0.3)]
        out_correct[judge_outrange]=np.nan
        logging.warning('Coordinates '+str(judge_outrange)+' are out of application range (E(B-V)>=0.3)')
        
        unvalid_index = unvalid_index[(out[unvalid_index]<0.3)]
        l_forfit = l[unvalid_index]
        b_forfit = b[unvalid_index]
        
        # Read Planck2014 dust temperature map
         # TAU:0,RADIANCE:3,TEMP:4,BETA:6
        planck_temp = hp.fitsfunc.read_map('data/HFI_CompMap_ThermalDustModel_2048_R1.20.fits',\
                                  field=4,nest=False,hdu=1)
        temp = planck_temp[indices[unvalid_index]]
        
        # Read Planck2014 beta map
        # TAU:0,RADIANCE:3,TEMP:4,BETA:6
        planck_beta = hp.fitsfunc.read_map('data/HFI_CompMap_ThermalDustModel_2048_R1.20.fits',\
                                  field=6,nest=False,hdu=1)
        beta = planck_beta[indices[unvalid_index]]

        ### using Planck2014-Tau fit function to get correction factor k ###
        
        # read the values of fitting parameter:
        func = pd.read_csv('data/kfit_relation/all_func.csv')
        Plank2014R_func = func['Planck2014-Tau'].values

        # get normed k from the fitting relation
        k_norm_fit = kfit_func([temp,beta],Plank2014R_func)

        # convert normed k_fit to k_fit based on k_median-EBV relation
        k2normk_file = pd.read_csv('data/k2normk/all_k2normk_pub.csv')
        k2normk_ebv = k2normk_file['Planck2014-Tau'].values
        k2normk_kmedian = k2normk_file['Planck2014-Tau_k'].values
        f=interpolate.interp1d(k2normk_ebv,k2normk_kmedian,kind='linear')
        k_median = f(out[unvalid_index])
        k_fit = k_norm_fit*k_median
        
        #get final corrected ebv
        out_correct[unvalid_index] = out[unvalid_index] * k_fit
        
        # check whether k_fit is reasonable value
        # if kfit is negative or too large. it cannot be used.
        judge_kfit = unvalid_index[((k_fit<0)|(k_fit>2))]        
        if len(judge_kfit)!=0:
            out_correct[judge_kfit]=np.nan
            logging.warning('Coordinates '+str(judge_kfit)+' cannot get a reasonable k_fit value.')

    return out,out_correct


def Planck2019Tau_reddening_correction(coordinate):
    '''
    This function is for correcting Planck2019-Tau extinction map. (Irfan et al. 2019)
    
    Args:
        1. coordinate('astropy.coordinates.SkyCoord'): The coordinates of target to correct extinction. 
        Only two options can be accepted : 'icrs' or 'galactic'.
        
    Returns:
        1. out('numpy.float64'): Original reddening value E(B-V) from the Planck2019-Tau map;
            
        2. out_correct('numpy.float64'): The corrected reddening value from maps. 
        
    Warnings:
        1. 'Coordinates [i] are out of application range (E(B-V)>=0.3)': This means this pixel is neither 
        within the footprint of LAMOST where can use k map correct the reddening directly, or within the 
        proper range of using the k_fit relationship.
        
        2. 'Coordinates [1] cannot get a reasonable k_fit value.': This means the final k_fit of this pxiel
        is out of [0, 2], which we think it is not a reasonable value for correction.
    '''
    
    # convert coordinates to galactic frame     
    if coordinate.frame.name == 'icrs':
        coordinate = coordinate.galactic
        l = coordinate.galactic.l.deg
        b = coordinate.galactic.b.deg
    else:
        l = coordinate.l.deg
        b = coordinate.b.deg
    
    # prepare the origianl reddening value from the Planck2019 map 
    planck_map = hp.fitsfunc.read_map('data/tau.fits') 
    out = np.full(len(l), np.nan, dtype='f4')
    out_correct = np.full(len(l), np.nan, dtype='f4')
    
   # convert Skycoord to indice in Planck2019-Tau HEALPix map 
    nside_planck = 2048
    theta = np.radians(90. - b)
    phi = np.radians(l)
    indices = hp.ang2pix(nside_planck, theta, phi)

    # 1.49*10**4 is the factor of the transformation between EBV and Tau
    out = planck_map[indices]*1.49*10**4
        
    #load the k map of Planck2019-Tau
    use_k_map = hp.fitsfunc.read_map\
    ('data/k_maps/k_Planck2019-Tau.fits', nest=False, dtype=None) 
    
    # convert Skycoord to indice k HEALPix map    
    nside = 64 
    npix = hp.nside2npix(nside)
    theta = np.radians(90. - b)
    phi = np.radians(l)
    indices = hp.ang2pix(nside, theta, phi)
    
    ### read the value from k map ###
    k = use_k_map[indices]

    # check whether k is in the footprint. If not, k will be obtained from fitted k relation
    valid = (k > -1e+30)
    if np.any(valid):
        out_correct[valid] = out[valid]*k[valid]
        
    un_valid = (k < -1e+30)
    if np.any(un_valid):
        unvalid_index = np.where(un_valid==True)[0]
        
        #check whether the coordinate is in the applicable range of EBV
        judge_outrange = unvalid_index[(out[unvalid_index]>=0.3)]
        out_correct[judge_outrange]=np.nan
        logging.warning('Coordinates '+str(judge_outrange)+' are out of application range (E(B-V)>=0.3)')
        
        unvalid_index = unvalid_index[(out[unvalid_index]<0.3)]
        l_forfit = l[unvalid_index]
        b_forfit = b[unvalid_index]
        
        # Read Planck2019 dust temperature map
        planck_temp = hp.fitsfunc.read_map('data/temp.fits')
        temp = planck_temp[indices[unvalid_index]]

        # Read Planck2019 beta map
        planck_beta = hp.fitsfunc.read_map('data/beta.fits')
        beta = planck_beta[indices[unvalid_index]]
        
        ### using Planck2019-Tau fit function to get correction factor k ###
        
        # read the values of fitting parameter:
        func = pd.read_csv('data/kfit_relation/all_func.csv')
        Plank2014R_func = func['Planck2019-Tau'].values
        
        # get normed k from the fitting relation
        k_norm_fit = kfit_func([temp,beta],Plank2014R_func)
        
        # convert normed k_fit to k_fit based on k_median-EBV relation
        k2normk_file = pd.read_csv('data/k2normk/all_k2normk_pub.csv')
        k2normk_ebv = k2normk_file['Planck2019-Tau'].values
        k2normk_kmedian = k2normk_file['Planck2019-Tau_k'].values
        f=interpolate.interp1d(k2normk_ebv,k2normk_kmedian,kind='linear')
        k_median = f(out[unvalid_index])
        k_fit = k_norm_fit*k_median
        
        # get final corrected ebv
        out_correct[unvalid_index] = out[unvalid_index] * k_fit
        
        # check whether k_fit is reasonable value
        # if kfit is negative or too large. it cannot be used.  
        judge_kfit = unvalid_index[((k_fit<0)|(k_fit>2))]      
        if len(judge_kfit)!=0:
            out_correct[judge_kfit]=np.nan
            logging.warning('Coordinates '+str(judge_kfit)+' cannot get a reasonable k_fit value.')

    return out,out_correct