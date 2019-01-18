# -*- coding: utf-8 -*-
import time
import pysm
import numpy as np
import healpy as hp
from pysm.nominal import models
from pysm.common import convert_units
from NPTFit import create_mask as cm

def Mask(sky_map):
    mask = np.empty_like(total);mask_i = cm.make_mask_total(nside,b_mask=True, b_deg_min = -10, b_deg_max = 10)
    for i in range(Nf):
        for j in range(3):
            mask[i][j]=mask_i
    map_masked = hp.ma(sky_map)
    map_masked.mask = mask
    return map_masked
    
def simulation(nside):
    c1_config = models("c1", nside)       
    s1_config = models("s1", nside)
    d1_config = models("d1", nside)
    f1_config = models("f1", nside)
    a1_config = models("a1", nside)
    sky_config = {
        'synchrotron' : s1_config,
        'dust' : d1_config,
        'freefree' : f1_config,
        'cmb' : c1_config,
        'ame' : a1_config
    }
    sky = pysm.Sky(sky_config) 
    
    global wn, cmb_MASKED, coefficients,mask,cmb
    nu = np.array([30., 70.,150., 353.]) 
    coefficients = convert_units("uK_RJ", "uK_CMB", nu)
     
    dust = sky.dust(nu);synchrotron = sky.synchrotron(nu);freefree = sky.freefree(nu)
    cmb = sky.cmb(nu); ame = sky.ame(nu)
    #np.save('../simulation/Nside = 1024/dust_map.npy',convert_unit(dust)); np.save('../simulation/Nside = 1024/synchro_map.npy',convert_unit(synchrotron));
    np.save('../simulation/Nside = 1024/4_fre/cmb_map.npy',convert_unit(cmb)); #np.save('../simulation/Nside = 1024/freefree_map.npy',convert_unit(freefree));
    #np.save('../simulation/Nside = 1024/ame_map.npy',convert_unit(ame))
    
    total = sky.signal()(nu)
    #total = (total - cmb)*0.7 + cmb         #decrease the foreground 
    np.save('../simulation/Nside = 1024/4_fre/total_map.npy',convert_unit(total))
    
def convert_unit(map):
    for i in range(0,Nf):
        map[i] = map[i]*coefficients[i]
    return map
    
def power_spectrum(a,R):
    '''calculate the power spectrum and re-organize as the form of R_Y(l) = A*R_S*A.T + R_n'''
    cl = []; Cl = []                                                                           #cl = (x00, x01 ,,, x32, x33)
    for i in range(Nf):                                                                        # x00[l=0,1,2,3,4,5,,,] x01 x02 x03 x(Nf)
        for j in range(Nf):                                                                    # x10 x11 x12 x13    # x20 x21 x22 x23                                                        
            c_pq = hp.anafast(a[i][0], a[j][0], lmax = L,gal_cut=R)        # x(Nf,0) x31 x32 x33,,,x(Nf,Nf)(l = 0,1,2,,,L),,, gal_cut=10
            cl.append(c_pq)                                                                                                                                                                                                                                                                                              
    
    for l in range(L+1):
        x = np.zeros(Nf*Nf)
        for i in range(0,Nf*Nf):
            x[i] = cl[i][l]     
        temp = x.reshape(Nf,Nf)
        Cl.append(temp)
    return Cl 

#### larissa 
#L = 1500
#Nf = 7
#total = np.load('/home/yao/Desktop/ABS/lylirsa/ABS/total_map.npy')
#total_ps = power_spectrum(total,0)
#np.save('/home/yao/Desktop/ABS/lylirsa/ABS/total_power_spectrum.npy', total_ps)
#exit()
####

if __name__ == '__main__':
    start = time.time()
    Nf = 4; L =2000; Q = 50;nside= 1024
    simulation(nside)
    ###false_number = len(mask[mask==False])
    ##fr = 1 #float(false_number)/(12*nside**2)
    total = np.load('../simulation/Nside = 1024/4_fre/total_map.npy')
    ###ame = np.load('../simulation/ame_map.npy')
    ##synchro = np.load('../simulation/Nside = 1024/synchro_map.npy')
    ##noise = np.load('../simulation/noise/noise_map.npy')    
    ##total = total + noise + synchro
    ##test_total_power = power_spectrum(total,10)
    ##np.save('../simulation/test/total_power_with_noise_synch*2',test_total_power)
    ##print 'hahaha'
    ##exit()
   
    #cmb = np.load('../simulation/Nside = 1024/shift/cmb_map.npy')
    ##total_masked = Mask(total)
    ##cmb_masked = Mask(cmb)
    ##CL = hp.anafast((cmb)[0][0],lmax = L,gal_cut=10)    #cmb_MASKED,gal_cut=10   
    ##np.save('../simulation/Nside = 1024/masked_outside/cmb_power_spectrum.npy', CL)    
    
    #white_noise
    for n in range(0,50):
        #white_noise = np.load('/home/yao/Desktop/EM-smica/ABS/simulation/Nside = 1024/noise_map_NO%s.npy'%(n+2))
        white_noise = np.ones_like(total); cl_noise = np.ones((Nf,L+1));  Nl = (0.066**2,0.063**2,0.015**2,0.068**2)#(0.066**2,0.065**2,0.063**2,0.028**2,0.015**2,0.023**2,0.068**2); 
        #Nl = np.array(Nl); Nl = 100*(Nl); print Nl
        Nl_real = []; Nl_matrix=[]
        for i in range(Nf):
            cl_noise[i] = Nl[i] 
            wn_j = hp.synfast(cl_noise[i],nside)
            for j in range(3):
                white_noise[i,j,:]= wn_j
                
            noise_power_spectrum = hp.anafast((white_noise)[i][0], lmax = L,gal_cut=10)
            Nl_real.append(noise_power_spectrum)
        Nl_real = np.array(Nl_real);Nl_T = Nl_real.T
        for i in range(L):
            Nl_matrix.append(np.diag((Nl_T[i])))
        #np.save('/home/yao/Desktop/EM-smica/ABS/simulation/noise/noise_map_NO%s.npy'%(n+2), white_noise)
        np.save('/home/yao/Desktop/EM-smica/ABS/simulation/Nside = 1024/4_fre/noise_power_spectrum_NO%s.npy'%(n+1), Nl_matrix)        
        
        #white_noise_masked = Mask(white_noise)
        #total_masked = Mask(total - cmb + white_noise)
    
        total_masked_power_spectrum = power_spectrum(total+white_noise,10)
        ###total_power_spectrum = power_spectrum(total + white_noise)         
        np.save('../simulation/Nside = 1024/4_fre/total_power_spectrum_NO%s.npy'%(n+1), total_masked_power_spectrum) 
        ###np.save('../simulation/total_power_spectrum_NO%s.npy'%(n+1), total_power_spectrum)
        print "You have arrive at step:%s"%n
    print "Simulation done!"      
    
    #r_n = np.mat(np.diag(Nl))   
    #R_n = []
    #for l in range(L+1):
        #R_n.append(r_n)
    #wn_power_spectrum = R_n
    #np.save('../simulation/noise_power_spectrum.npy', wn_power_spectrum)

    end = time.time()
    print "Running time is %s minutes"%((end - start)/60)
