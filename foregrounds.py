import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from NPTFit import create_mask as cm

total = np.load('/home/yao/Desktop/EM-smica/ABS/simulation/Nside = 1024/total_map.npy')
Nf=7
def Mask(sky_map):
    mask = np.empty_like(total);mask_i = cm.make_mask_total(1024,b_mask=True, b_deg_min = -10, b_deg_max = 10)
#     mask_i  = cm.make_mask_total(1024,band_mask=True, band_mask_range = 10)
    for i in range(Nf):
        for j in range(3):
            mask[i][j]=mask_i
    map_masked = hp.ma(sky_map)
    map_masked.mask = mask
    return map_masked


dust = np.load('/home/yao/Desktop/EM-smica/ABS/simulation/Nside = 1024/dust_map.npy')
synch =  np.load('/home/yao/Desktop/EM-smica/ABS/simulation/Nside = 1024/synchro_map.npy')
freefree =  np.load('/home/yao/Desktop/EM-smica/ABS/simulation/Nside = 1024/freefree_map.npy')
noisepower =  np.load('/home/yao/Desktop/EM-smica/ABS/simulation/Nside = 1024/masked_outside/noise_power_spectrum_real.npy')
cmb = np.load('/home/yao/Desktop/EM-smica/ABS/simulation/Nside = 1024/cmb_map.npy')
ame = np.load('/home/yao/Desktop/EM-smica/ABS/simulation/Nside = 1024/ame_map.npy')

dust = Mask(dust);synch=Mask(synch);freefree=Mask(freefree);cmb=Mask(cmb);ame=Mask(ame)
names = [30, 44, 70, 95, 150, 217, 353]

L = 2000;R=0
noise = [] 
for i in range(7):
    noise_i = []
    for l in range(L):
        noise_i.append(noisepower[l][i,i])
    noise.append(noise_i)
    
for n in range(7):
    dust_power = hp.anafast(dust[n][0],lmax=L,gal_cut=R)
    synch_power = hp.anafast(synch[n][0],lmax=L,gal_cut=R)
    freefree_power = hp.anafast(freefree[n][0],lmax=L,gal_cut=R)
    cmb_power = hp.anafast(cmb[n][0],lmax=L,gal_cut=R)
    ame_power = hp.anafast(ame[n][0],lmax=L,gal_cut=R)
    ell = np.arange(len(dust_power))
    
    plt.loglog(ell[1:],(ell[1:]*(ell[1:]+1)*cmb_power[1:]/2/np.pi),label = 'cmb',lw=2)
    plt.loglog(ell[1:],(ell[1:]*(ell[1:]+1)*dust_power[1:]/2/np.pi),label = 'dust',lw=2)
    plt.loglog(ell[1:],(ell[1:]*(ell[1:]+1)*synch_power[1:]/2/np.pi),label = 'synch',lw=2)
    plt.loglog(ell[1:],(ell[1:]*(ell[1:]+1)*freefree_power[1:]/2/np.pi),label = 'freefree',lw=2)
    plt.loglog(ell[1:],(ell[1:]*(ell[1:]+1)*ame_power[1:]/2/np.pi),label = 'ame',lw=2)
    nll = np.arange(len(noise[n]))
    plt.loglog(nll[1:],(nll[1:]*(nll[1:]+1)*noise[n][1:]/2/np.pi),label = 'noise',lw=2)
    plt.xticks([1e4,1e3,1e2,1e1,1e0],[r'$10000$',r'$1000$',r'$100$',r'$10$',r'$1$'])
    #plt.plot(ell,ell*(ell+1)*zll/2/np.pi,label = '150')
    plt.xlabel(r'$\ell$', fontsize=16)
    plt.ylabel(r'$\ell(\ell+1)C_\ell$/2$\pi$ [$\mu$k$^2$]',fontsize=16)

    plt.legend(loc = 'lower right',frameon = False)
    plt.title('%s GHz'%names[n])
    plt.savefig("/home/yao/Desktop/EM-smica/ABS/simulation/Nside = 1024/%s_GHz_fore_mask_outside.pdf"%names[n],format='pdf')
    plt.close()
    #plt.show()
    print n