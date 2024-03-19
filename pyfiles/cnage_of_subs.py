import fitsio
import numpy as np
import matplotlib.pyplot as plt

allstars = fitsio.read('./catalogs/horta_stars_dr17.fits')

# n20all = np.where(allstars['SUB_NAME']=='Sequoia(N20)')[0]
# print(len(n20all))

# print(len(allstars))
# print(allstars['FE_H_ERR'])
# print(max(allstars['LOGG']),min(allstars['LOGG']))
# print(max(allstars['FE_H']),min(allstars['FE_H']))

# names = ['Aleph','Arjuna','GES','HelmiStream','Heracles','Nyx','Sagittarius','Sequoia(M19)','Sequoia(N20)','Thamnos','Iitoi']
# fehlim = np.where(allstars['FE_H']>-1.5)[0]
# fehlimstars = allstars[fehlim]
# print(len(fehlim))
# for name in names:
#     numstars = np.where(fehlimstars['SUB_NAME']==name)[0]
#     print(name,len(numstars))

usable = np.where((allstars['FE_H']>-1.2) & (allstars['FE_H']<.3) & (allstars['LOGG']>2.14) & (allstars['LOGG']<3.1))[0]
# print(len(usable))

gdstars = allstars[usable]
#print(gdstars)

class RGiant(object):

    # a class of red giant stars that are aware of useful parameters + catalog data

    def __init__(self,row=None,cfe=None,nfe=None,feh=None,logg=None,teff=None,
                 gaiaID=None,sdssID=None,subname=None,hortapath=None):

        if hortapath is None:
            self.hortapath = './catalogs/horta_stars_dr17.fits'
        else:
            pass
        if row is not None:
            self.subname = row['SUB_NAME']
            self.cfe = row['C_FE']
            self.nfe = row['N_FE']
            self.feh = row['FE_H']
            self.logg = row['LOGG']
            self.teff = row['TEFF']

            self.gaiaEDR3ID = row['GAIAEDR3_SOURCE_ID']
            self.sdssID = row['APOGEE_ID']
            self.row = row
        else:
            print('No star data')

        self.cnr = None
        self.age = None

    def get_cnr(self):
        self.cnr = self.cfe - self.nfe
        return self.cnr
        
    def get_age(self):
        m = 2.412
        b = 10.204
        self.age = m * self.cnr + b
        return self.age
        
stars = []
for row in gdstars:
    star = RGiant(row=row)
    star.get_cnr()
    star.get_age()
    # print(star.subname,star._cnr,star._age)
    stars.append(star)
    # print(star.subname)

gdstars = []
for star in stars:
    if star.age is not np.isnan(star.age):
        gdstars.append(star)

subnames = []
ages = []
fehs = []

for star in gdstars:
    ages.append(star.age)
    # print(10**star.age/1e9)
    subnames.append(star.subname)
    fehs.append(star.feh)

ages = np.array(ages)
fehs = np.array(fehs)
subnames = np.array(subnames)
unisubnames = np.unique(subnames)
np.save('ages',ages)
np.save('subnames',subnames)
np.save('fehs',fehs)
print(max(ages),min(ages))

# counts,edges,bars = plt.hist(ages)

# plt.savefig('./plots/hist_test.pdf')

# for name in unisubnames:
#     mems = np.where(subnames == name)[0]
    
#     streamage = ages[mems]
#     streamfeh = fehs[mems]
#     counts,edges,bars = plt.hist(streamage)
#     plt.savefig('./plots/hist_test_'+str(name)+'.pdf')
#     plt.clf()

#     plt.scatter(streamage,streamfeh)
#     plt.savefig('./plots/test_'+str(name)+'.pdf')
#     plt.clf()


# for i in ages:
#     gyr = (10**i)/1e9
#     print(gyr,i)


arjuna = []
arjages = []
arjfeh = []
ges = []
gesages = []
gesfeh = []
helmi = []
helages = []
helfeh = []
heracles = []
herages = []
herfeh = []
nyx = []
nyxages = []
nyxfeh = []
sag = []
sagages = []
sagfeh = []
m19 = []
m19ages = []
m19feh = []
n20 = []
n20ages = []
n20feh = []
thamnos = []
thamages = []
thamfeh= []

arjcnr = []
gescnr = []
helcnr = []
hercnr = []
nyxcnr = []
sagcnr = []
m19cnr = []
n20cnr = []
thamcnr = []

for star in gdstars:
    if star.subname == unisubnames[1]:
        arjuna.append(star)
        arjages.append(star.age)
        arjfeh.append(star.feh)
        arjcnr.append(star.cnr)
    elif star.subname == unisubnames[2]:
        ges.append(star)
        gesages.append(star.age)
        gesfeh.append(star.feh)
        gescnr.append(star.cnr)
    elif star.subname == unisubnames[3]:
        helmi.append(star)
        helages.append(star.age)
        helfeh.append(star.feh)
        helcnr.append(star.cnr)
    elif star.subname == unisubnames[4]:
        heracles.append(star)
        herages.append(star.age)
        herfeh.append(star.feh)
        hercnr.append(star.cnr)
    elif star.subname == unisubnames[6]:
        nyx.append(star)
        nyxages.append(star.age)
        nyxfeh.append(star.feh)
        nyxcnr.append(star.cnr)
    elif star.subname == unisubnames[7]:
        sag.append(star)
        sagages.append(star.age)
        sagfeh.append(star.feh)
        sagcnr.append(star.cnr)
    elif star.subname == unisubnames[9]:
        m19.append(star)
        m19ages.append(star.age)
        m19feh.append(star.feh)
        m19cnr.append(star.cnr)
    elif star.subname == unisubnames[10]:
        n20.append(star)
        n20ages.append(star.age)
        n20feh.append(star.feh)
        n20cnr.append(star.cnr)
    elif star.subname == unisubnames[11]:
        thamnos.append(star)
        thamages.append(star.age)
        thamfeh.append(star.feh)
        thamcnr.append(star.cnr)

# print(len(arjuna),len(arjages))
# print(len(ges),len(gesages))
# print(len(helmi),len(helages))
# print(len(heracles))
# print(len(nyx))
# print(len(sag))
# print(len(m19))
# print(len(n20))
# print(len(thamnos),len(thamages))
ages_list = [(arjages,arjcnr),(gesages,gescnr),(helages,helcnr),(herages,hercnr),(nyxages,nyxcnr),(sagages,sagcnr),(m19ages,m19cnr),(n20ages,n20cnr),(thamages,thamcnr)]

edges = np.linspace(7,12,21)
# print(edges)
# counts,bins = np.histogram(arjages,edges,[7.,12.])
# plt.stairs(counts,bins,color='forestgreen',fill=True)
# plt.savefig('./plots/stairs_test.pdf')


clist,blist = [],[]
for item in ages_list:
    counts,bins = np.histogram(item[0],edges,[7.,12.])
    clist.append(counts)
    blist.append(bins)

fig, ((ax1,ax2,ax3),(ax4,ax5,ax6),(ax7,ax8,ax9)) = plt.subplots(3,3,sharex='row')
fig.set_size_inches(21,16)
# fig, ((ax1,ax2,ax3)) = plt.subplots(1,3)
# fig.suptitle('Distribution of Ages by Sub-structure',fontsize=30)
# fig.supxlabel('log(age(yr))')
# fig.supylabel('Count')

# ax1.hist(arjages,color='forestgreen')
ax1.stairs(clist[0],blist[0],color='forestgreen',fill=True)
ax1.set_title('Arjuna',fontsize=24)
ax1.tick_params(axis='both',labelsize=18)
ax1.vlines(10.24,-2,300,linestyles='dashed',color='black',label='log(age) of Universe')
ax1.set_ylim(0,8.5)
ax1.legend(loc='upper left',fontsize=18)

# ax2.hist(gesages,color='cornflowerblue')
ax2.stairs(clist[1],blist[1],color='cornflowerblue',fill=True)
ax2.set_title('GES',fontsize=24)
ax2.tick_params(axis='both',labelsize=18)
ax2.vlines(10.24,-2,300,linestyles='dashed',color='black')
ax2.set_ylim(0,175)

# ax3.hist(helages,color='firebrick')
ax3.stairs(clist[2],blist[2],color='firebrick',fill=True)
ax3.set_title('HelmiStream',fontsize=24)
ax3.tick_params(axis='both',labelsize=18)
ax3.vlines(10.24,-2,300,linestyles='dashed',color='black')
ax3.set_ylim(0,8.5)

# ax4.hist(herages,color='darkorchid')
ax4.stairs(clist[3],blist[3],color='darkorchid',fill=True)
ax4.set_title('Heracles',fontsize=24)
ax4.tick_params(axis='both',labelsize=18)
ax4.vlines(10.24,-2,300,linestyles='dashed',color='black')
ax4.set_ylim(0,27.5)

# ax5.hist(nyxages,color='pink')
ax5.stairs(clist[4],blist[4],color='pink',fill = True)
ax5.set_title('Nyx',fontsize=24)
ax5.tick_params(axis='both',labelsize=18)
ax5.vlines(10.24,-2,300,linestyles='dashed',color='black')
ax5.set_ylim(0,44.5)

# ax6.hist(sagages,color='goldenrod')
ax6.stairs(clist[5],blist[5],color='goldenrod',fill=True)
ax6.set_title('Sagittarius',fontsize=24)
ax6.tick_params(axis='both',labelsize=18)
ax6.vlines(10.24,-2,300,linestyles='dashed',color='black')
ax6.set_ylim(0,8.5) 

# ax7.hist(m19ages,color='peru')
ax7.stairs(clist[6],blist[6],color='peru',fill=True)
ax7.set_title('Sequoia(M19)',fontsize=24)
ax7.tick_params(axis='both',labelsize=18)
ax7.vlines(10.24,-2,300,linestyles='dashed',color='black')
ax7.set_ylim(0,9.5)

# ax8.hist(n20ages,color='teal')
ax8.stairs(clist[7],blist[7],color='teal',fill=True)
ax8.set_title('Sequoia(N20)',fontsize=24)
ax8.tick_params(axis='both',labelsize=18)
ax8.vlines(10.24,-2,300,linestyles='dashed',color='black')
ax8.set_ylim(0,8.5)

# ax9.hist(thamages,color='slategrey')
ax9.stairs(clist[8],blist[8],color='slategrey',fill=True)
ax9.set_title('Thamnos',fontsize=24)
ax9.tick_params(axis='both',labelsize=18)
ax9.vlines(10.24,-2,300,linestyles='dashed',color='black')
ax9.set_ylim(0,12.5)

fig.savefig('./plots/comb/comb_hist.pdf',format='pdf',bbox_inches='tight',pad_inches=.05)


# def plot_stuff(axis,x,y,c,barcolor,out,y_err = None,color='k'):
#     outlier = np.logical_and(~np.isnan(x), ~np.isnan(y))
#     for i in y[outlier]:
#         if np.abs(i) >= out:
#             print(f"{c.name} has a star with a delta value of {i}")
    
#     axis.hlines(0,-300,7500,colors = 'k',alpha = 0.8,linestyle = 'dashed')
#     axis.set_xlim(min(x[outlier])-.01,max(x[outlier])+.01)
#     #axis.xaxis.set_major_locator(plt.MaxNLocator(5))
#     # axis.errorbar(x[outlier], y[outlier], yerr=y_err[outlier],ecolor='grey',fmt='none', capsize=5, zorder=0)
#     axis.errorbar(-0.2,-0.4, yerr = y_err,xerr = y_err,c='tab:grey',fmt='.', capsize=5, zorder=0,markersize=1)
#     axis.scatter(x[outlier],y[outlier],c=barcolor,cmap = 'viridis', s=35,vmin=1,vmax=5)
#     image = axis.scatter(x[outlier],y[outlier],c=barcolor,cmap = 'viridis', s=35,vmin=1,vmax=5) #change vmin,vmax if doing alphas
#     axis.text(min(x[outlier]),40,"{:}".format(c.name), fontsize = 'xx-large')  #.13 -> .1 for alphas

#     return image
plt.clf()
# plt.scatter(arjages,arjfeh,c='forestgreen',s=25)

# plt.scatter(gesages,gesfeh,c='cornflowerblue',s=25)

# plt.scatter(helages,helfeh,c='firebrick',s=25)

# plt.scatter(herages,herfeh,c='darkorchid',s=25)

# plt.scatter(nyxages,nyxfeh,c='pink',s=25)


# plt.savefig('./plots/comb/age_feh.pdf')

fig1, ((ax1,ax2,ax3),(ax4,ax5,ax6),(ax7,ax8,ax9)) = plt.subplots(3,3,sharex=True,sharey=True)
fig1.set_size_inches(24,24)

ax1.scatter(arjages,arjfeh,c='forestgreen',s=30)
ax1.xaxis.set_tick_params(which='both',labelbottom=True)
ax1.tick_params(axis='both',labelsize=21)
ax1.set_title('Arjuna',fontsize=30)
# ax1.errorbar(7.5,-1,yerr=.015,xerr=.05,fmt='none')


ax2.scatter(gesages,gesfeh,c='cornflowerblue',s=30)
ax2.xaxis.set_tick_params(which='both',labelbottom=True)
ax2.yaxis.set_tick_params(which='both',labelleft=True)
ax2.tick_params(axis='both',labelsize=21)
ax2.set_title('GES',fontsize=30)

ax3.scatter(helages,helfeh,c='firebrick',s=30)
ax3.xaxis.set_tick_params(which='both',labelbottom=True)
ax3.yaxis.set_tick_params(which='both',labelleft=True)
ax3.tick_params(axis='both',labelsize=21)
ax3.set_title('HelmiStream',fontsize=30)

ax4.scatter(herages,herfeh,c='darkorchid',s=30)
ax4.xaxis.set_tick_params(which='both',labelbottom=True)
ax4.tick_params(axis='both',labelsize=21)
ax4.set_title('Heracles',fontsize=30)

ax5.scatter(nyxages,nyxfeh,c='pink',s=30)
ax5.xaxis.set_tick_params(which='both',labelbottom=True)
ax5.yaxis.set_tick_params(which='both',labelleft=True)
ax5.tick_params(axis='both',labelsize=21)
ax5.set_title('Nyx',fontsize=30)

ax6.scatter(sagages,sagfeh,c='goldenrod',s=30)
ax6.xaxis.set_tick_params(which='both',labelbottom=True)
ax6.yaxis.set_tick_params(which='both',labelleft=True)
ax6.tick_params(axis='both',labelsize=21)
ax6.set_title('Sagittarius',fontsize=30)

ax7.scatter(m19ages,m19feh,c='peru',s=30)
ax7.tick_params(axis='both',labelsize=21)
ax7.set_title('Sequoia(M19)',fontsize=30)

ax8.scatter(n20ages,n20feh,c='teal',s=30)
ax8.yaxis.set_tick_params(which='both',labelleft=True)
ax8.tick_params(axis='both',labelsize=21)
ax8.set_title('Sequoia(N20)',fontsize=30)

ax9.scatter(thamages,thamfeh,c='slategrey',s=30)
ax9.yaxis.set_tick_params(which='both',labelleft=True)
ax9.tick_params(axis='both',labelsize=21)
ax9.set_title('Thamnos',fontsize=30)



fig1.savefig('./plots/comb/age_feh.pdf',format='pdf',bbox_inches='tight',pad_inches=.05)



names = ['Arjuna','GES','HelmiStream','Heracles','Nyx','Sagittarius','Sequoia(M19)','Sequoia(N20)','Thamnos']

for stream,name in zip(ages_list,names):
    stream_arr = np.array(stream)
    nans = np.isnan(stream_arr[0])
    avg = np.average(stream_arr[0][~nans])
    med = np.median(stream_arr[0][~nans])
    std = np.std(stream_arr[0][~nans])
    nanscnr = np.isnan(stream_arr[1])
    cnravg = np.average(stream_arr[1][~nanscnr])
    medcnr = np.median(stream_arr[1][~nanscnr])
    stdcnr = np.std(stream_arr[1][~nanscnr])
    print(name+' Age:')
    print('~~~~~'*4)
    print('Median: ',np.round(med,2))
    print('Average: ',np.round(avg,2))
    print('Standard Deviation: ',np.round(std,2))
    print('~~~~~'*4)
    print(name+' [C \ N]:')
    print('Median: ',np.round(medcnr,2))
    print('Average: ',np.round(cnravg,2))
    print('Standard Deviation: ',np.round(stdcnr,2))
    print('~~~~~'*4)

