import fitsio

stars = fitsio.read('./catalogs/horta_stars_dr17.fits')
stars['GAIA_ID']

stars.dtype.names