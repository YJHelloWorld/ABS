import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

def bin_l(cl):
    bin_averages = []
    for l in range(L):
        cl[l] = l*(l+1)/2/np.pi*(cl[l])    
    for q in range(Q):
        bin_averages.append(sum(cl[q*L/Q:((q+1)*L/Q)]/(L/Q)))
    return bin_averages
        
L = 2000; Q=50; Nf=7; E_cut= 0.5; D_B = []; SamNum = 14; Evals = np.ones((Q,Nf)); Evals_all = [];Delta = 51*10

for n in range(SamNum):
    D_B_n = []; 
    loc = '/home/yao/Desktop/EM-smica/ABS/simulation/Nside = 1024/null_test/b_40'
    #D =  np.load('/home/yao/Desktop/EM-smica/ABS/simulation/Nside = 1024/unmasked/noise_power_spectrum_besides_diag.npy')
    #D = np.load('/home/yao/Desktop/EM-smica/ABS/simulation/test/total_power_with_noise_synch*2.npy')
    D = np.load('%s/total_power_spectrum_NO%s.npy'%(loc,n+1))
    noise_ps = np.load('%s/noise_power_spectrum_real.npy'%(loc))
    #CMB = np.load('/home/yao/Desktop/EM-smica/ABS/simulation/Nside = 1024/null_test/cmb_power_spectrum.npy')  
    D = bin_l(D); noise_ps = bin_l(noise_ps)  #*fr ; CMB = bin_l(CMB)
    noise = np.zeros_like(noise_ps);f = []
    noise = np.load('%s/noise_power_spectrum_RMS_2000.npy'%(loc))
    noise = bin_l(noise)
    noise = np.array(noise)
    #print noise.max()
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
                D[l][i,j] = D[l][i,j] /np.sqrt(noise[l][i]*noise[l][j]) + Delta*f[l][i]*f[l][j] 
                
    for l in range(Q): 
        e_vals,E = LA.eig(D[l])
        #print e_vals
        Evals[l,:] = e_vals 
        
    ##begin of eigen_value test            
    #fig3 = plt.figure(3)
    #evals = Evals.T; evals_nega = np.ones_like(evals)
    #for i in range(Nf):
        #for j in range(Q):
            #if evals[i,j]<= 0:
                #evals_nega[i,j] = abs(evals[i,j])
                #evals[i,j] = None
            #else:
                #evals_nega[i,j] = None
    #ell = np.ones(Q)
    #for q in range(Q):
        #ell[q] = (2*q+1)*L/Q/2    
    #for i in range(Nf):  
        #x = np.arange(len(evals[i]))
        #plt.scatter(ell,(evals[i][:]),color = 'g') #label = 'positive eigenvalues'
        #plt.scatter(ell,evals_nega[i][:],color = 'r') #,label = 'negative eigenvalues'
        #plt.axhline(0.5,color = 'k')
        #plt.yscale('log')
    #plt.xlabel(r'$\ell$');plt.ylabel('Eigenvalues')
    #plt.savefig('/home/yao/Desktop/EM-smica/ABS/simulation/Nside = 1024/null_test/Eigenvalues_%s'%n,format = 'pdf')
    #plt.close()
#exit()
    ##end of eigen_value test
        
        for i in range(Nf):
            E[:,i]=E[:,i]/LA.norm(E[:,i])**2  
       
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
#frame1=fig1.add_axes((.1,.3,.8,.6))
#ax = plt.gca()
#majors = np.linspace(200,1000,5)
#ax.yaxis.set_major_locator(ticker.FixedLocator(majors))
for n in range(SamNum):
    ell = np.ones(Q)
    for q in range(Q):
        ell[q] = (2*q+1)*L/Q/2
    if n ==0:
        plt.plot(ell,D_B[n],'g-^',label = 'recovered CMB') #ell*(ell+1)/2/np.pi*
        
    else:
        plt.plot(ell,D_B[n],'g-^') #ell*(ell+1)/2/np.pi*
plt.axhline(0,color = 'k')
#plt.plot(ell,ell*(ell+1)*CMB[0:Q]/(2*np.pi),'r-x',label = 'real CMB')
plt.xlabel(r'$\ell$'); plt.ylabel(r'$\ell(\ell+1)C_\ell$/2$\pi$ [$\mu$k$^2$]')
#plt.ylim(0,1100)
#plt.ylim(100,5000)
plt.legend()
#begin error bar
#frame1.set_xticklabels([])
#frame2=fig1.add_axes((.1,.1,.8,.2))
#ax = plt.gca()
#majors = np.linspace(-0.03,0.01,7)
#ax.yaxis.set_major_locator(ticker.FixedLocator(majors))
fig2 = plt.figure(2)
D_B_whole = np.array(D_B)
DB = D_B_whole.T
sigma_A = 0; delta_A = 0
Ell = np.ones(Q)
SE = []; Mean = []
for i in range(0,Q):
    Ell[i] = (2*i+1)*L/Q/2 
#     DB[i] = ell[i]*(ell[i]+1)*DB[i]
    se = np.std(DB[i]) # /CMB[i]  *Ell*(Ell+1)/2/np.pi
    SE.append(se)
    #sigma_A = sigma_A + (np.average(DB[i])/se)**2
    #delta_A = delta_A + np.average(DB[i])*CMB[i]/se**2
    mean = np.average(DB[i]) #/CMB[i]-1
    Mean.append(mean)
    plt.errorbar(Ell[i],mean,yerr=se,fmt='-o')
    plt.axhline(0,color = 'k')
data = np.column_stack((Mean,SE))
np.savetxt('/home/yao/Desktop/EM-smica/ABS/codes/null_test_b_40_data',data,fmt="%s")
#print "delta_A = ", delta_A/sigma_A-1
#sigma_A = np.sqrt(1/sigma_A)
#print "sigma_A = ", sigma_A
plt.ylabel(r'$\Delta$ $\ell(\ell+1)C_\ell$/2$\pi$  ')
#plt.savefig('/home/yao/Desktop/EM-smica/ABS/simulation/Nside = 1024/null_test/recovered CMB without first point',format = 'pdf')
plt.legend()
#end error bar

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
##plt.savefig('/home/yao/Desktop/EM-smica/ABS/simulation/Nside = 1024/null_test/Eigenvalues',format = 'pdf')

plt.legend()
plt.show()
