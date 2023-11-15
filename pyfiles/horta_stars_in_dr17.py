## PYTHON lIBS ##
from astropy.io import fits
from astropy.table import Table, Column
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
import numpy.lib.recfunctions as rec
import os

# read in the dr17 file
with fits.open('/scratch/jotto/sdssIV/allStar-dr17-synspec_rev1.fits') as hdul:
    # hdul.info()
    dr17 = hdul[1].data


# CLEAN THE DATA UP! Removing bad bits and -9999 values for APOGEE. Grabbing only cluster members from OCCAM.
# TWO BITWISE FLAGS FOR BAD DATA             
badbits = 2**23        # aspcapstar flag - Chemistry
suspectbits = 2**16    # star flag - Stellar parameters
cfe_bad = 2**21        # aspcap flag - bad [C/Fe] values
nfe_bad = 2**22        # aspcap flag - bad [N/Fe] values
snr_bad = 2**27        # aspcap flag - BAD S/N flag (S/N<50)

gd = (np.bitwise_and(dr17['ASPCAPFLAG'], badbits) == 0) &\
(np.bitwise_and(dr17['STARFLAG'], suspectbits) == 0) &\
(np.bitwise_and(dr17['ASPCAPFLAG'],cfe_bad ) == 0) &\
(np.bitwise_and(dr17['ASPCAPFLAG'], nfe_bad) == 0) &\
(np.bitwise_and(dr17['ASPCAPFLAG'], snr_bad) == 0) & (dr17['SNREV']>75)
teff_logg_check = np.logical_and(dr17["TEFF"] > 0, dr17["LOGG"] > -10)
teff_logg_feh_check = np.logical_and(dr17["FE_H"] > -6, teff_logg_check)
C_N_check = np.logical_and(dr17["N_FE"]>-2000, dr17["CI_FE"] > -2000) &\
(dr17["C_FE_ERR"]<=0.1) & (dr17["N_FE_ERR"]<=0.1) & (dr17["C_FE"] > -2000)
evol_cut = (dr17["LOGG"]<3.3) & (dr17['LOGG']>1)
#star_param_check = np.logical_and(teff_logg_feh_check, CI_N_check)
indices1 = np.where(np.logical_and(gd,teff_logg_feh_check), evol_cut, C_N_check)
#cleaned allstar date. no -9999 values
good_apo = dr17[indices1]

# print(len(dr17), len(good_apo))
# print(f"APOGEE - CLEAN APOGEE = {len(dr17)-len(good_apo)} removed")

with fits.open('./catalogs/halo_subs_Horta.fits') as hdul:
    # hdul.info()
    horta_halos = hdul[1].data

# print(np.unique(horta_halos['name']))

apo_ind = []
halo_ind = []

for i, halostar in enumerate(horta_halos['APOGEE_ID']):
    star_id = np.where(good_apo['APOGEE_ID']==halostar)[0]
    if len(star_id) > 0:
        apo_ind.append(star_id[0])
        halo_ind.append(i)

print(len(apo_ind),len(halo_ind))

substruc_names = np.core.records.fromarrays([horta_halos[halo_ind]['name']],names='SUB_NAME')

horta_subhalo = rec.merge_arrays([dr17[apo_ind],substruc_names],flatten=True)
print(len(horta_subhalo))

hdu = fits.BinTableHDU(data=horta_subhalo)
print(hdu.columns)
# hdulist = fits.HDUList([hdu])
# hdu.writeto('./catalogs/horta_stars_dr17.fits')

#print(horta_subhalo)

