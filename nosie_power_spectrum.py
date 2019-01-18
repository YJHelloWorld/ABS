import numpy as np
import healpy as hp

fr = 1;nu = [30, 95, 150, 217, 353];Nf = 5
ps = (0.066**2*fr,0.028**2*fr,0.015**2*fr,0.023**2*fr,0.068**2*fr)
cl = np.ones((5,1501));WN=[];nl = []
for i in range(4,5):
    for j in range(200):
        wn = hp.synfast(cl[i]*ps[i],512)
        WN.append(wn)
        nl_j = hp.anafast(wn,lmax=1500)
        nl.append(nl_j) 
        print j
    np.save('/home/yao/Desktop/EM-smica/ABS/simulation/noise/WNmap_%sGHz_512_100.npy'%nu[i],WN)
    np.save('/home/yao/Desktop/EM-smica/ABS/simulation/noise/Nl_%sGHz_512_200.npy'%nu[i],nl)