import fitsio
import numpy as np

allstars = fitsio.read('./catalogs/horta_stars_dr17.fits')
print(len(allstars))

usable = np.where((allstars['FE_H']>-1.2) & (allstars['FE_H']<.3) & (allstars['LOGG']>2.14) & (allstars['LOGG']<3.1))[0]

gdstars = allstars[usable]
print(len(gdstars))
#print(gdstars)

class GalSubStruc(object):

    def __init__(self, name=None, memstarIDs=None, memstarsdata=None, cfes=None ,nfes=None, fehs=None,
                 hortapath = None):

        self.name = name
        self.memstarIDs = memstarIDs
        self.memstarsdata = memstarsdata
        self.cfes = cfes
        self.nfes = nfes
        self.fehs = fehs
        self.hortapath = hortapath

        self.cnrs = None
        self.ages = None

        if self.memstarIDs is None:
            print("No stars for this structure. That's not possible.")

        if self.memstarsdata is None:
            print("No data for the structure")

        if self.cfes is None:
            print("No carbon")

        if self.nfes is None:
            print("No nitrogen")
        
        if self.fehs is None:
            print("No iron")

        if self.hortapath is None:
            print("Missing path to stars")

    def get_cnrs(self):

        if not (len(self.cfes) == len(self.nfes)):
            print("carbon and nitrgen have different numbers of stars")
        gd = np.isnan(self.cfes) | np.isnan(self.nfes)
            
        self.cnrs = self.cfes[~gd] - self.nfes[~gd] 
        return self.cnrs
      
    def get_ages(self):

        m = 2.41
        b = 10.20
        self.ages = m * self.cnrs + b
        return self.ages

def makeStructure(name,cat):
    
    stars = np.where(cat['SUB_NAME']==name)[0]
    strucStars = cat[stars]

    apIDs = strucStars['APOGEE_ID']
    gaiaIDS = strucStars['GAIAEDR3_SOURCE_ID']
    IDs = np.column_stack((apIDs,gaiaIDS))

    carbon = strucStars['C_FE']
    nitrogen = strucStars['N_FE']
    iron = strucStars['FE_H']

    path = './catalogs/horta_stars_dr17.fits'

    return GalSubStruc(name=name,memstarIDs=IDs,memstarsdata=strucStars,cfes=carbon,
                       nfes=nitrogen,fehs=iron,hortapath=path)



ges = makeStructure('GES',gdstars)
ges.get_cnrs()
ges.get_ages()
# print(ges.cnrs)

# print(ges.ages)
# print('~~~'*4)
# print(ges.cnrs)

ages = (ges.ages**10)/1e9
# print(ages)
print(max(ages))
print(min(ages))
print(np.median(ages))
print(np.sum(ages)/len(ages))