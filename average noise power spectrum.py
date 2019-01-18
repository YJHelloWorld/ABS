import numpy as np
#import matplotlib.pyplot as plt
import healpy as hp

SamNum = 14; Nf = 7;L=2000;Q=2000
Nl = np.zeros((L,Nf,Nf));                   #  When calculating the rms values of noise power spectrum, these should be commented out.
nl=[];Nl_std=[]

def bin_l(cl):
    bin_averages = []
    for l in range(L):
        cl[l] = l*(l+1)/2/np.pi*(cl[l])    
    for q in range(Q):
        bin_averages.append(sum(cl[q*L/Q:((q+1)*L/Q)]/(L/Q)))
    return bin_averages

loc = '/home/yao/Desktop/EM-smica/ABS/simulation/Nside = 1024/null_test/b_20'

for n in range(SamNum):
    n_diag=[]
    nl_i = np.load('%s/noise_power_spectrum_NO%s.npy'%(loc, n+1))
    ###nl_i = bin_l(nl_i) ##
    #for i in range(Q): ##
        #n_diag.append(np.diag(nl_i[i])) ##   When calculating the mean values of noise power spectrum, i.e. the theoretical values, should be commented out.
    Nl = Nl + nl_i                           #
    
    #nl.append(n_diag) ##
#for i in range(Q): ##
    #nl_s=[];nl_std=np.empty(Nf) ##
    #for n in range(SamNum): ##
        #nl_s.append(nl[n][i]) ##
    #tran = np.array(nl_s).T   ##
    #for f in range(Nf): ## 
        #nl_std[f] = np.std(tran[f]) ##
    #Nl_std.append(nl_std) ##
    
#np.save('%s/noise_power_spectrum_RMS_2000.npy'%(loc),Nl_std) ##
 
Nl = Nl/SamNum                                 #
np.save('%s/noise_power_spectrum_real.npy'%loc,Nl)  #

