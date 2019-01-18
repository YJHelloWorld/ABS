import numpy as np
import healpy as hp

total = np.load('../simulation/Nside = 1024/noise_map.npy')
for n in range(20):
    #white_noise = np.load('/home/yao/Desktop/EM-smica/ABS/simulation/Nside = 1024/noise_map_NO%s.npy'%(n+2))
    white_noise = np.ones_like(total); cl_noise = np.ones((Nf,L+1)); Nl = (0.066**2,0.065**2,0.063**2,0.028**2,0.015**2,0.023**2,0.068**2); 
    Nl_real = []; Nl_matrix=[]
    for i in range(Nf):
        cl_noise[i] = Nl[i] 
        wn_j = hp.synfast(cl_noise[i],nside)
        for j in range(3):
            white_noise[i,j,:]= wn_j
    
        noise_power_spectrum = hp.anafast(white_noise[i][0], lmax = L,gal_cut=0)
        Nl_real.append(noise_power_spectrum)
    Nl_real = np.array(Nl_real);Nl_T = Nl_real.T
    for i in range(L):
        Nl_matrix.append(np.diag((Nl_T[i])))
    #np.save('/home/yao/Desktop/EM-smica/ABS/simulation/noise/noise_map_NO%s.npy'%(n+2), white_noise)
    np.save('/home/yao/Desktop/EM-smica/ABS/simulation/Nside = 1024/noise_power_spectrum_NO%s.npy'%(n+1), Nl_matrix)        
    
    ##white_noise_masked = Mask(white_noise)
    #total_masked_power_spectrum = power_spectrum(total + white_noise,0)
    ##total_power_spectrum = power_spectrum(total + white_noise)         
    #np.save('../simulation/total_power_spectrum_NO%s.npy'%(n+2), total_masked_power_spectrum) 
    ##np.save('../simulation/total_power_spectrum_NO%s.npy'%(n+1), total_power_spectrum)
    print "You have arrive at step:%s"%n
print "Simulation done!" 