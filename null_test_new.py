#import sys 
#sys.path.append('/home/yao/Desktop/EM-smica/smica-codes/programs') 
#from smica_simulation_new import power_spectrum
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from NPTFit import create_mask as cm
import healpy as hp
#mask = cm.make_mask_total(512,b_mask=True, b_deg_min = -10, b_deg_max = 10) #leave the middle unmasked
mask = cm.make_mask_total(512,band_mask=True, band_mask_range = 10) #mask the middle 
z = np.ones(len(mask))
fr = 1 -np.dot(mask,z)/len(mask) #masked
print fr

def bin_l(cl):
    bin_averages = []
    for l in range(L):
        cl[l] = l*(l+1)/2/np.pi*(cl[l])    
    for q in range(Q):
        bin_averages.append(sum(cl[q*L/Q:((q+1)*L/Q)]/(L/Q)))
    return bin_averages
        
L = 1000; Q=50; Nf=7; E_cut= 0.5; D_B = []; SamNum = 50; Evals = np.ones((Q,Nf)); Evals_all = [];Delta = 18*10#.0037*10

for n in range(SamNum):
    D_B_n = []; 
    #D =  np.load('/home/yao/Desktop/EM-smica/ABS/simulation/Nside = 1024/unmasked/noise_power_spectrum_besides_diag.npy')
    #D = np.load('/home/yao/Desktop/EM-smica/ABS/simulation/test/total_power_with_noise_synch*2.npy')
    D = np.load('/home/yao/Desktop/EM-smica/ABS/simulation/Nside = 1024/null_test/512/total_power_spectrum_0.1cmb_NO%s.npy'%(n+1))
    noise_ps = np.load('/home/yao/Desktop/EM-smica/ABS/simulation/Nside = 1024/null_test/512/noise_power_spectrum_real_0.1cmb.npy')
    #cmb_map = np.load('../simulation/Nside = 1024/cmb_map_512.npy')
    #CMB = hp.anafast((0.1*cmb_map)[0][0],lmax = L,gal_cut=10)    #cmb_MASKED,gal_cut=10   
    CMB = np.load('/home/yao/Desktop/EM-smica/ABS/simulation/Nside = 1024/null_test/512/cmb_power_spectrum.npy')      
    D = bin_l(D);  noise_ps = bin_l(noise_ps)  ;CMB = bin_l(CMB);
    noise = np.zeros_like(noise_ps);f = []
    noise = np.load('/home/yao/Desktop/EM-smica/ABS/simulation/Nside = 1024/null_test/512/noise_power_spectrum_RMS_0.1cmb.npy')
    noise = bin_l(noise)
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
    ## begin eigenvalues test
    #for n in range(SamNum):
        #ell = np.ones(Q)
        #for q in range(Q):
            #ell[q] = (2*q+1)*L/Q/2    
    #if n ==0:
        #fig3 = plt.figure(3)
        #evals = Evals.T; evals_nega = np.ones_like(evals)
        #for i in range(Nf):
            #for j in range(Q):
                #if evals[i,j]<= 0:
                    #evals_nega[i,j] = abs(evals[i,j])
                    #evals[i,j] = None
                #else:
                    #evals_nega[i,j] = None
        #for i in range(Nf):  
            #x = np.arange(len(evals[i]))
            #plt.scatter(ell,(evals[i][:]),color = 'g') #label = 'positive eigenvalues'
            #plt.scatter(ell,evals_nega[i][:],color = 'r') #,label = 'negative eigenvalues'
            #plt.axhline(0.5,color = 'k')
            #plt.yscale('log')
        #plt.xlabel(r'$\ell$');plt.ylabel('Eigenvalues')
        #plt.savefig('/home/yao/Desktop/EM-smica/ABS/results/eigval_Nf=7_lbin=30_masked_realization%s.png'%n,format = 'png')
        #plt.show()
        ##end of test
        
        D_B_l = 0
        for i in range(Nf):
            if e_vals[i]>=E_cut:
                G_i = np.dot(f[l],E[:,i])
                D_B_l = D_B_l + (G_i**2/e_vals[i]) 
        D_B_l = 1/ D_B_l - Delta
        D_B_n.append(D_B_l)
       
    #Evals_all.append(Evals)
    D_B.append(D_B_n)
#np.save('/home/yao/Desktop/EM-smica/ABS/simulation/masked_outside/D_B_masked_outside.npy',D_B)

fig1 = plt.figure(1)  
frame1=fig1.add_axes((.1,.3,.8,.6))
#ax = plt.gca()
#majors = np.linspace(200,1000,5)
#ax.yaxis.set_major_locator(ticker.FixedLocator(majors))
for n in range(SamNum):
    ell = np.ones(Q)
    for q in range(Q):
        ell[q] = (2*q+1)*L/Q/2
    if n ==0:
        plt.plot(ell,D_B[n],'g-^',label = 'recovered CMB') #ell*(ell+1)/(2*np.pi)
    else:
        plt.plot(ell,D_B[n],'g-^')  #ell*(ell+1)/(2*np.pi)
plt.plot(ell,CMB[0:Q],'r-x',label = 'real CMB') #ell*(ell+1)/(2*np.pi)
plt.xlabel(r'$\ell$'); plt.ylabel(r'$\ell(\ell+1)C_\ell/2\pi$ [$\mu$k$^2$]')
#plt.ylim(20,5000)
#plt.ylim(100,5000)
plt.legend()
#begin error bar
frame1.set_xticklabels([])
frame2=fig1.add_axes((.1,.1,.8,.2))
#ax = plt.gca()
#majors = np.linspace(-0.03,0.01,7)
#ax.yaxis.set_major_locator(ticker.FixedLocator(majors))
D_B_whole = np.array(D_B)
DB = D_B_whole.T
sigma_A = 0; delta_A = 0
Ell = np.ones(Q)
SE = []
for i in range(0,Q):
    Ell[i] = (2*i+1)*L/Q/2 
#     DB[i] = ell[i]*(ell[i]+1)*DB[i]
    se = np.std((DB[i]-CMB[i]))/CMB[i]
    SE.append(se)
       
    
    #sigma_A = sigma_A + (np.average(DB[i])/se)**2
    #delta_A = delta_A + np.average(DB[i])*CMB[i]/se**2
    mean = np.average((DB[i]-CMB[i]))/CMB[i]
    plt.errorbar(Ell[i],mean,yerr=se,fmt='-o')
    plt.axhline(0,color = 'k')
print np.mean(SE)
np.savetxt('null_test_error',SE,fmt="%s")
#print "delta_A = ", delta_A/sigma_A-1
#sigma_A = np.sqrt(1/sigma_A)
#print "sigma_A = ", sigma_A
plt.ylabel(r'$\Delta C_\ell / C_\ell^{real} - 1$ [%]')
#plt.yticks([-0.020,-0.015,-0.010, -0.005,0,0.005,0.010],[-2,-1.5,-1.0,-0.5,0,0.5,1.0])
#plt.savefig('/home/yao/Desktop/EM-smica/ABS/simulation/Nside = 1024/null_test/recovered CMB without first point',format = 'pdf')
plt.legend()
#end error bar

#fig2 = plt.figure(2)
#Error = [];s=2
#for i in range(s,Q):
    #error = (D_B[0][i]-CMB[i])/CMB[i]*100
    #Error.append(error)
#plt.plot(ell[s:Q],Error)
#plt.xlabel(r'$\ell$');plt.ylabel('Fractional error /%')
##plt.savefig('/home/yao/Desktop/EM-smica/ABS/simulation/masked_outside/masked_outside_error_lbin=30',format = 'pdf')
#####np.save('/home/yao/Desktop/EM-smica/ABS/results/Errors_1',Error)

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
#plt.savefig('/home/yao/Desktop/EM-smica/ABS/simulation/Nside = 1024/null_test/Eigenvalues',format = 'pdf')

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

