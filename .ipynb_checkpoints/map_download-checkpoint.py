import requests
import sys

def url_download(link, file_name, save_dir):
    
    with open(save_dir+file_name, "wb") as f:
        print("Downloading %s" % file_name)
        response = requests.get(link, stream=True, allow_redirects=True)
        total_length = response.headers.get('content-length')

        if total_length is None: # no content length header
            f.write(response.content)
        else:
            dl = 0
            total_length = int(total_length)
            for data in response.iter_content(chunk_size=4096):
                dl += len(data)
                f.write(data)
                done = int(50 * dl / total_length)
                sys.stdout.write("\r[%s%s]" % ('=' * done, ' ' * (50-done)) )    
                sys.stdout.flush()

                
def SFD_query():
    
    save_dir = 'data/'
    
    # download dust maps
    url = 'https://dataverse.harvard.edu/api/access/datafile/:persistentId?persistentId=doi:10.7910/DVN/EWCNL5/JOV0CM'
    file_name = 'SFD_dust_4096_ngp.fits'
    url_download(url, file_name, save_dir)
    
    url = 'https://dataverse.harvard.edu/api/access/datafile/:persistentId?persistentId=doi:10.7910/DVN/EWCNL5/XTLXX5'
    file_name = 'SFD_dust_4096_sgp.fits'
    url_download(url, file_name, save_dir)
    
    # download temp maps
    url = 'https://dataverse.harvard.edu/api/access/datafile/:persistentId?persistentId=doi:10.7910/DVN/EWCNL5/FVLE7F'
    file_name = 'SFD_temp_ngp.fits'
    url_download(url, file_name, save_dir)
    
    url = 'https://dataverse.harvard.edu/api/access/datafile/:persistentId?persistentId=doi:10.7910/DVN/EWCNL5/D52DVN'
    file_name = 'SFD_temp_sgp.fits'
    url_download(url, file_name, save_dir)
    
    
def Planck2014_query():
    
    save_dir = 'data/'
    
    # download Planck2014 maps
    url = 'https://irsa.ipac.caltech.edu/data/Planck/release_1/all-sky-maps/maps/HFI_CompMap_ThermalDustModel_2048_R1.20.fits'
    file_name = 'HFI_CompMap_ThermalDustModel_2048_R1.20.fits'
    url_download(url, file_name, save_dir)
    
    
def Planck2019_query():
    
    save_dir = 'data/'
    
    # download Planck2019 maps
    url = 'http://cdsarc.u-strasbg.fr/ftp/cats/J/A+A/623/A21/fits/beta.fits'
    file_name = 'beta.fits'
    url_download(url, file_name, save_dir)
    
    url = 'http://cdsarc.u-strasbg.fr/ftp/cats/J/A+A/623/A21/fits/temp.fits'
    file_name = 'temp.fits'
    url_download(url, file_name, save_dir)

    url = 'http://cdsarc.u-strasbg.fr/ftp/cats/J/A+A/623/A21/fits/tau.fits'
    file_name = 'tau.fits'
    url_download(url, file_name, save_dir)