
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from NPTFit import create_mask as cm
#mask = cm.make_mask_total(512,b_mask=True, b_deg_min = -10, b_deg_max = 10) #leave the middle unmasked
mask = cm.make_mask_total(512,band_mask=True, band_mask_range = 10) #mask the middle 
z = np.ones(len(mask))
fr = 1 -np.dot(mask,z)/len(mask) #masked

def bin_l(cl):
    global q_cut,l_cut
    bin_averages = []
    for q in range(Q):
        if q*bin1 < 900:
            bin_averages.append(sum(cl[q*bin1:((q+1)*bin1)]/(bin1)))
        elif (q-1)*bin1 < 900:
            q_cut = q; l_cut = q_cut*bin1 
        elif l_cut + (q-1 - q_cut)*bin2 < L:
            print l_cut + (q -1- q_cut)*bin2
            bin_averages.append(sum(cl[(l_cut+(q-1-q_cut)*bin2):(l_cut+(q-q_cut)*bin2)]/(bin2)))
    return bin_averages
        
L = 1500; Q=50; bin1= 30; bin2 = 120; Nf=7; E_cut=1; D_B = []; SamNum = 1;  Delta = 0.00434117706763*20
q_cut = 0; l_cut = 0

for n in range(SamNum):
    D_B_n = []; 
    D = np.load('/home/yao/Desktop/EM-smica/ABS/simulation/total_power_spectrum_masked_NO%s.npy'%(n+1))
    noise_ps = np.load('/home/yao/Desktop/EM-smica/ABS/simulation/noise_power_spectrum.npy')
    CMB = np.load('/home/yao/Desktop/EM-smica/ABS/simulation/cmb_power_spectrum_masked.npy')    
    D = bin_l(D); CMB = bin_l(CMB); noise_ps = bin_l(noise_ps*fr)  #*fr
    Q = len(D)
    Evals = np.ones((Q,Nf))
    noise = np.zeros_like(noise_ps);f = []
    for i in range(Q):
        f_q = np.ones(Nf)
        if i <30:
            bins = bin1
        else:
            bins = bin2
        for j in range(Nf):
            noise[i][j,j] = noise_ps[i][j,j]*np.sqrt(2.0/((2*i+1)*fr)*bins) # *fr
            f_q[j] = f_q[j]/np.sqrt(noise[i][j,j])
        f.append(f_q)
          
    for l in range(Q):
        D[l] = D[l] - noise_ps[l]
        for i in range(Nf): 
            for j in range(Nf):
                D[l][i,j] = D[l][i,j]/np.sqrt(noise[l][i,i]*noise[l][j,j]) + Delta*f[l][i]*f[l][j] 
                
    for l in range(Q): 
        e_vals,E = LA.eig(D[l])
        if n == 0:
            Evals[l,:] = e_vals
        for i in range(Nf):
            E[:,i]=E[:,i]/LA.norm(E[:,i])**2  
        
        D_B_l = 0
        for i in range(Nf):
            if e_vals[i]>=E_cut:
                G_i = np.dot(f[l],E[:,i])
                D_B_l = D_B_l + (G_i**2/e_vals[i]) 
        D_B_l = 1/ D_B_l - Delta
        D_B_n.append(D_B_l)
    D_B.append(D_B_n)
#np.save('/home/yao/Desktop/EM-smica/ABS/simulation/test/D_B_0.5.npy',D_B)

fig1 = plt.figure(1)  
#frame1=fig1.add_axes((.1,.3,.8,.6))
print q_cut
for n in range(SamNum):
    ell = np.ones(Q)
    for q in range(Q):
        if q < q_cut:
            ell[q] = (2*q+1)*bin1/2
        else:
            ell[q] = l_cut + (2*(q-q_cut)+1)*bin2/2
            print len(ell),Q
    plt.plot(ell,ell*(ell+1)*D_B[n]/(2*np.pi),'g-^') #,label = 'our CMB'    
plt.plot(ell,ell*(ell+1)*CMB[0:Q]/(2*np.pi),'r-x',label = 'real CMB')
plt.savefig('/home/yao/Desktop/EM-smica/ABS/simulation/test/cmb_masked_with_diff_lbin_revised',format = 'eps')

#frame1.set_xticklabels([])
##begin error bar
#frame2=fig1.add_axes((.1,.1,.8,.2))
#ax = plt.gca()
#majors = np.linspace(-0.008,0.004,7)
#ax.yaxis.set_major_locator(ticker.FixedLocator(majors))
#D_B_whole = np.load('/home/yao/Desktop/EM-smica/ABS/results/D_B_0.5.npy');DB = D_B_whole.T
#sigma_A = 0; delta_A = 0
#Ell = np.ones(Q)
#for i in range(1,Q):
    #Ell[i] = (2*i+1)*L/Q/2 
##     DB[i] = ell[i]*(ell[i]+1)*DB[i]
    #se = np.std(DB[i]/CMB[i]) #
    #sigma_A = sigma_A + (DB[i][0]/se)**2
    #delta_A = delta_A + DB[i][0]*CMB[i]/se**2
    #mean = np.average(DB[i])/CMB[i]-1
    #plt.errorbar(Ell[i],mean,yerr=se,fmt='-o')
    #plt.axhline(0,color = 'k')
#print "sigma_A = ", sigma_A
#print "delta_A = ", delta_A/sigma_A-1
#print delta_A
#plt.savefig('/home/yao/Desktop/EM-smica/ABS/results/cmb_unmasked_with_error',format = 'eps')
#plt.legend()
##end error bar

fig2 = plt.figure(2)
Error = [];s=1
for i in range(s,Q):
    error = (D_B[0][i]-CMB[i])/CMB[i]*100
    Error.append(error)
plt.plot(ell[s:Q],Error)
plt.savefig('/home/yao/Desktop/EM-smica/ABS/simulation/test/Errors_lbin_revised',format = 'pdf')

fig3 = plt.figure(3)
evals = Evals.T
for i in range(Nf):
    for j in range(Q):
        if evals[i,j]<= 0:
            evals[i,j] = None 
for i in range(Nf):  
    x = np.arange(len(evals[i]))
    plt.scatter(ell,(evals[i][:]))
    plt.axhline(0.5,color = 'k')
    plt.yscale('log')
plt.savefig('/home/yao/Desktop/EM-smica/ABS/simulation/test/Eignvalues_lbin_revised',format = 'eps')

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
