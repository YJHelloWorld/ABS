
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt

db_masked_0_1 = np.load('/home/yao/Desktop/EM-smica/ABS/results/Errors_0.1_masked.npy')
db_0_1 = np.load('/home/yao/Desktop/EM-smica/ABS/results/Errors_0.1.npy')
db_masked_1 = np.load('/home/yao/Desktop/EM-smica/ABS/results/Errors_1_masked.npy')
db_1 = np.load('/home/yao/Desktop/EM-smica/ABS/results/Errors_1.npy')

L=1500;Q=50
ell = np.ones(Q)
for q in range(Q):
    ell[q] = (2*q+1)*L/Q/2

plt.plot(ell[1:Q],db_masked_0_1,'-o',label = 'Masked_0.1')
plt.plot(ell[1:Q],db_0_1,'-^',label = 'Unasked_0.1')
plt.plot(ell[1:Q],db_masked_1,'-s',label = 'Masked_1')
plt.plot(ell[1:Q],db_1,'-v',label = 'Unasked_1')
plt.ylabel('%')
plt.legend()
plt.savefig("/home/yao/Desktop/EM-smica/ABS/results/Errors",format='svg')
plt.show()