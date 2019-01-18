#import sys 
#sys.path.append('/home/yao/Desktop/EM-smica/smica-codes/programs') 
#from smica_simulation_new import power_spectrum
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from NPTFit import create_mask as cm
#mask = cm.make_mask_total(512,b_mask=True, b_deg_min = -10, b_deg_max = 10) #leave the middle unmasked
mask = cm.make_mask_total(512,band_mask=True, band_mask_range = 10) #mask the middle 
z = np.ones(len(mask))
fr = 1 -np.dot(mask,z)/len(mask) #masked
print fr

def bin_l(cl):
    bin_averages = []
    for q in range(Q):
        bin_averages.append(sum(cl[q*L/Q:((q+1)*L/Q)]/(L/Q)))
    return bin_averages
        
L = 2000; Q=50; Nf=7; E_cut= 0.5; D_B = []; SamNum = 50; Evals = np.ones((Q,Nf)); Evals_all = [];Delta = 0#.0014*10

for n in range(1):
    D_B_n = []; 
    D =  np.load('/home/yao/Desktop/EM-smica/ABS/simulation/Nside = 1024/unmasked/noise_power_spectrum_besides_diag.npy')   
    #D = np.load('/home/yao/Desktop/EM-smica/ABS/simulation/test/total_power_with_noise_synch*2.npy')   
    #D = np.load('/home/yao/Desktop/EM-smica/ABS/simulation/total_power_spectrum_NO%s.npy'%(n+1))
    noise_ps = np.load('/home/yao/Desktop/EM-smica/ABS/simulation/Nside = 1024/unmasked/noise_power_spectrum_real_test.npy')
    CMB = np.load('/home/yao/Desktop/EM-smica/ABS/simulation/cmb_power_spectrum.npy')  
    D = bin_l(D); CMB = bin_l(CMB); noise_ps = bin_l(noise_ps)  #*fr
    noise = np.zeros_like(noise_ps);f = []
    noise = np.load('/home/yao/Desktop/EM-smica/ABS/simulation/Nside = 1024/unmasked/noise_power_spectrum_RMS.npy')
    for i in range(Q):
        f_q = np.ones(Nf)
        for j in range(Nf):
            #noise[i][j,j] = noise_ps[i][j,j]*np.sqrt(2.0/(((2*i+1)*fr)*L/Q)) # *fr
            f_q[j] = f_q[j]/np.sqrt(noise[i][j]) #,j
        f.append(f_q)    
    for l in range(Q):
        D[l] = D[l] - noise_ps[l]
        for i in range(Nf): 
            for j in range(Nf):
                D[l][i,j] = D[l][i,j]/np.sqrt(noise[l][i]*noise[l][j]) + Delta*f[l][i]*f[l][j] 
        #if l ==48:
            #print np.diag(D[l])
                
    for l in range(Q): 
        e_vals,E = LA.eig(D[l])
        #if n ==0:
        Evals[l,:] = e_vals        
        
        for i in range(Nf):
            E[:,i]=E[:,i]/LA.norm(E[:,i])**2  
   
for n in range(SamNum):
    ell = np.ones(Q)
    for q in range(Q):
        ell[q] = (2*q+1)*L/Q/2 

fig3 = plt.figure(3)
evals = Evals.T; evals_nega = np.ones_like(evals)
for i in range(Nf):
    for j in range(Q):
        if evals[i,j]<= 0:
            evals_nega[i,j] = abs(evals[i,j])
            evals[i,j] = None
        else:
            evals_nega[i,j] = None
for i in range(Nf):  
    x = np.arange(len(evals[i]))
    plt.scatter(ell,(evals[i][:]),color = 'g') #label = 'positive eigenvalues'
    plt.scatter(ell,evals_nega[i][:],color = 'r') #,label = 'negative eigenvalues'
    plt.axhline(0.5,color = 'k')
    plt.yscale('log')
plt.xlabel(r'$\ell$');plt.ylabel('Eigenvalues')
#plt.savefig('/home/yao/Desktop/EM-smica/ABS/simulation/masked_outside/masked_outside_eignvalues_lbin=30',format = 'pdf')

#fig4 = plt.figure(4)
#z = []
#for i in range(50):
    #z.append(Evals[i].min())
#plt.plot(ell,z)
#plt.title('The minimum eigenvalues of each l_bin')

#cl = np.ones(1500)*0.028**2
#zll = np.ones(1500)*0.015**2
#ell = np.arange(1500)
#plt.plot(ell,ell*(ell+1)*cl/2/np.pi,label = '95')
#plt.plot(ell,ell*(ell+1)*zll/2/np.pi,label = '150')
plt.legend()
plt.show()

