import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
from scipy.integrate import quad


# Best-fit Planck values arXiv:1807.06209
#H0 = 67.66
#OM = 0.3111
#OL = 0.6889
#OK = 0
c = 299792.458


f = open("SN_data_sorted", "r")
l = f.readlines()

z = []
mu = []
sigmamu = []

for x in l:
	z.append(float(x.split(' ')[0]))
	mu.append(float(x.split(' ')[1]))
	sigmamu.append(float(x.split(' ')[2]))

z = np.asarray(z)
mu = np.asarray(mu)
sigmamu = np.asarray(sigmamu)

'''
dL = c/70 * (1+z) * np.sin(np.log(1+z))

with open('SN_dL_data.txt', 'a') as f:
	for zz, dd in zip(z, dL):
		f.write("%f %f\n" % (zz, dd))

'''

sigma_ext = 0.01
sigma_sample = 0.15

sigma = np.sqrt(sigmamu**2 + sigma_ext**2 + sigma_sample**2)


Harray = np.linspace(66.25, 66.3, num=101)
OMarray = np.linspace(0.279, 0.280, num=101)


def chisqfunc_timesphere(h): # d_L = c/H0 (1+z) sin(ln(1+z))
	a = h
	dL = c/a * (1+z) * np.sin(np.log(1+z))
	model = -5 + 5*np.log10(dL*1000000)
	chisq = np.sum( ((mu - model)/sigma)**2 )
	return chisq
	
def chisqfunc_rhct(h): # d_L = c/H0 (1+z) ln(1+z)
	a = h
	dL = c/a * (1+z) * np.log(1+z)
	model = -5 + 5*np.log10(dL*1000000)
	chisq = np.sum( ((mu - model)/sigma)**2 )
	return chisq
	
def E_z(z, OM, OK, OL):
    return 1/np.sqrt((1 + z)**3 * OM + (1 + z)**2 * OK + OL)
	
def chisqfunc_lcdm(arr):
	H0 = arr[0]
	Omega_m0 = arr[1]
	Omega_l0 = arr[2]
	Omega_k0 = 0
	
	temp = []
	for i in range(len(z)):
		temp.append((1+z[i]) * (c/H0) * quad(E_z, 0, z[i], args=(Omega_m0, Omega_k0, Omega_l0))[0])
	dL = np.asarray(temp)	
	model = -5 + 5*np.log10(dL*1000000)
	chisq = np.sum( ((mu - model)/sigma)**2 )
	return chisq
	
res_ts = []
res_rhct = []
res_lcdm = np.zeros((len(Harray), len(OMarray)))

for i in range(len(Harray)):
	res_ts.append(chisqfunc_timesphere(Harray[i]))
	res_rhct.append(chisqfunc_rhct(Harray[i]))
	
	#for j in range(len(OMarray)):
		#res_lcdm[i][j] = chisqfunc_lcdm((Harray[i], OMarray[j], 1-OMarray[j]))

res_ts = np.asarray(res_ts)
res_rhct = np.asarray(res_rhct)


#for h, cs_th, cs_rh in zip(Harray, res_ts, res_rhct):
    #print("%f %f %f" % (h, cs_th, cs_rh))
    

(inda, indb) = np.unravel_index(res_lcdm.argmin(), res_lcdm.shape)
#print(Harray[inda], OMarray[indb], res_lcdm[inda][indb])
    
print(Harray[np.argmin(res_ts)], min(res_ts))
print(Harray[np.argmin(res_rhct)], min(res_rhct))


'''
fig, ax = plt.subplots(1)
plt.errorbar(z, mu, yerr=sigma, fmt='o', zorder=0)
plt.plot(z, -5 + 5*np.log10(  (c/Harray[np.argmin(res_ts)] * (1+z) * np.sin(np.log(1+z))) *1000000), label='Timesphere', zorder=5)
plt.plot(z, -5 + 5*np.log10(  (c/Harray[np.argmin(res_rhct)] * (1+z) * np.log(1+z)) *1000000), label='RH=cT', zorder=10)
plt.show()
'''
