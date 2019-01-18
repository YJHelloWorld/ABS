#import sys 
#sys.path.append('/home/yao/Desktop/EM-smica/smica-codes/programs') 
#from smica_simulation_new import power_spectrum
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker
from NPTFit import create_mask as cm

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
        
L = 2000; Q=50; Nf=9; E_cut= 0.5; D_B = []; SamNum = 30; Evals = np.ones((Q,Nf)); Evals_all = [];Delta = 0#51*50
loc = '/home/yao/Desktop/EM-smica/ABS/simulation/Nside = 1024/9_fre'
for n in range(SamNum):
    D_B_n = []; 
    #D =  np.load('/home/yao/Desktop/EM-smica/ABS/simulation/Nside = 1024/unmasked/noise_power_spectrum_besides_diag.npy')
    #D = np.load('/home/yao/Desktop/EM-smica/ABS/simulation/test/total_power_with_noise_synch*2.npy')
    D = np.load('%s/total_power_spectrum_NO%s.npy'%(loc,n+1))
    noise_ps = np.load('%s/noise_power_spectrum_real.npy'%(loc))
    CMB = np.load('%s/cmb_power_spectrum.npy'%(loc))  
    D = bin_l(D); CMB = bin_l(CMB); noise_ps = bin_l(noise_ps)  #*fr
    noise = np.zeros_like(noise_ps);f = []
    noise = np.load('%s/noise_power_spectrum_RMS_2000.npy'%(loc))
    noise = bin_l(noise)
    #print noise
    for i in range(Q):
        f_q = np.ones(Nf)
        for j in range(Nf):
            #noise[i][j,j] = noise_ps[i][j,j]*np.sqrt(2.0/(((2*i+1)*fr)*L/Q)) # *fr
            f_q[j] = f_q[j]/np.sqrt(noise[i][j]) #,j
            #print noise[i][j]
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
        #if n == 49:
            #if l ==2:
                #print e_vals
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
mpl.rcParams['image.cmap'] = 'jet'
fig1 = plt.figure(1)
plt.title('Galatic- plane cut for 9 channels',fontsize = 20)
plt.axis('off')
frame1=fig1.add_axes((.1,.4,.8,.5))
#frame1.axes.spines['top'].set_linewidth(1.5) 


#ax = plt.gca()
#majors = np.linspace(200,1000,5)
#ax.yaxis.set_major_locator(ticker.FixedLocator(majors))
#import matplotlib
#matplotlib.rcParams['image.cmap'] = 'jet'

#Sig = []
#for i in range (Q):
    #sig = 0
    #for n in range(7):
        #sig_i = 1/noise[i][n]
        #sig = sig + sig_i
    #Sig.append(1/sig)
#plt.plot(np.arange(50), Sig)
##plt.ylim(0.01,1.4)

for n in range(SamNum):
    ell = np.ones(Q)
    for q in range(Q):
        ell[q] = (2*q+1)*L/Q/2
    if n ==0:
        plt.plot(ell,D_B[n],'g-^',label = 'recovered CMB') #ell*(ell+1)/(2*np.pi)
    else:
        plt.plot(ell,D_B[n],'g-^')  #ell*(ell+1)/(2*np.pi)
plt.plot(ell,CMB[0:Q],'r-x',label = 'real CMB') #ell*(ell+1)/(2*np.pi)
plt.xlabel(r'$\ell$'); plt.ylabel(r'$\mathcal{D}_{\ell}$ [$\mu$k$^2$]')
plt.ylim(20,5000)
#plt.ylim(20,6000)
plt.legend()
#begin error bar
frame1.set_xticklabels([])
frame2=fig1.add_axes((.1,.1,.8,.3))
#ax = plt.gca()
#majors = np.linspace(-0.03,0.01,7)
#ax.yaxis.set_major_locator(ticker.FixedLocator(majors))
D_B_whole = np.array(D_B)
DB = D_B_whole.T
sigma_A = 0; delta_A = 0
Ell = np.ones(Q)
SE = []; SE_rela =[]; Mean=[] ; Mean_rela = []; Mean_DB=[]

#loca = '/home/yao/Desktop/EM-smica/ABS/simulation/Nside = 1024'
#Cl_in = np.load('%s/masked/cmb_power_spectrum.npy'%loca)
#Cl_out = np.load('%s/masked_outside/cmb_power_spectrum.npy'%loca)
#cl_in = np.array(bin_l(Cl_in))
#cl_out = np.array(bin_l(Cl_out))
#re_cl_in = np.array(np.loadtxt('%s/masked/recovered_cl'%loca))
#re_cl_out = np.array(np.loadtxt('%s/masked_outside/recovered_cl'%loca))
#error_in = np.array(np.loadtxt('%s/masked/error_bar'%loca))
#error_out = np.array(np.loadtxt('%s/masked_outside/error_bar'%loca))

for i in range(0,Q):
    Ell[i] = (2*i+1)*L/Q/2 
#     DB[i] = ell[i]*(ell[i]+1)*DB[i]
    se = np.std((DB[i]-CMB[i]))/CMB[i]
    #print se
    SE.append(se)
    se_rela = np.std(DB[i]-CMB[i])#/CMB[i]
    SE_rela.append(se_rela)
    #sigma_A = sigma_A + (np.average(DB[i])/se)**2
    #delta_A = delta_A + np.average(DB[i])*CMB[i]/se**2
    mean = np.average((DB[i])-CMB[i])/CMB[i]
    mean_DB = np.average(DB[i])
    Mean_DB.append(mean_DB)
    mean_rela = np.average(DB[i]-CMB[i])#/CMB[i]
    Mean.append(mean)
    Mean_rela.append(mean_rela)
    plt.errorbar(Ell[i],100*mean,yerr=100*se,fmt='-o',capthick = 1.0)
    
plt.axhline(0,color = 'k')
#np.savetxt('/home/yao/Desktop/EM-smica/ABS/simulation/Nside = 1024/shift/shift_50.txt', Mean_DB)
#np.savetxt('%s/recovered_cl'%loc, (Mean_DB))
print Mean[0]
#np.savetxt('%s/error_bar'%loc,SE)

###plt.errorbar(ell,(((re_cl_in - cl_in)+ (re_cl_out - cl_out))/(cl_in+cl_out)) ,
             ###yerr = (error_in*0.826 + error_out*0.174),fmt='y-d',ecolor='y',capthick = 1.5) 

#plt.ylim( -0.06,0.085)
plt.ylim( -2.5,2.5)
z = np.column_stack((Mean_rela, SE_rela,Mean,SE))
np.savetxt('mean_and_error_rela_for_9fre', z,fmt="%s")

#print "delta_A = ", delta_A/sigma_A-1
#sigma_A = np.sqrt(1/sigma_A)
#print "sigma_A = ", sigma_A
plt.ylabel(r'$\Delta \mathcal{D}_\ell$/$\mathcal{D}^{real}_{\ell} $ [%]');plt.xlabel(r'$\ell$')
#plt.yticks([-0.060,-0.04,-0.020,0,0.020,0.04,0.06,0.08],[-6.0,-4.0,-2.0,0.0,2.0,4.0,6.0,8.0])
plt.yticks([-02,-1.5,-1,-0.5,-0.0,0.5,1,1.5,2],[-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0])
#plt.savefig('/home/yao/Desktop/EM-smica/ABS/simulation/Nside = 1024/null_test/recovered CMB without first point',format = 'pdf')
plt.legend()
#plt.savefig('/home/yao/Desktop/EM-smica/ABS/paper_figures/9_fre',format = 'pdf')
plt.show()
#end error bar
exit()

fig2 = plt.figure(2)
nll = np.arange(len(SE))
#plt.loglog(nll,SE)
#np.savetxt("stddev_CMB",SE)                
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
    plt.scatter(ell,evals_nega[i][:],color = 'r',marker='^') #,label = 'negative eigenvalues'
    plt.axhline(0.5,color = 'k')
    plt.yscale('log')
plt.xlabel(r'$\ell$');plt.ylabel(r'$\tilde{\lambda}_{\mu}$')
plt.xlim(-100,2100)
#plt.ylim(1e-3,1e7)
#plt.savefig('/home/yao/Desktop/EM-smica/ABS/simulation/Nside = 1024/2018-4-4/syn*2_masked_Eigenvalues.pdf',format = 'pdf')

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
