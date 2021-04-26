import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
import scipy.integrate as integrate

# Best-fit Planck values arXiv:1807.06209
#H0 = 67.66
#OM = 0.3111
#OL = 0.6889
#OK = 0
c = 299792.458


f = open("Hz_1807.07323.txt", "r")
l = f.readlines()[1:]

z = []
Hz = []
sigma = []

for x in l:
	z.append(float(x.split(' ')[0]))
	Hz.append(float(x.split(' ')[1]))
	sigma.append(float(x.split(' ')[2]))

z = np.asarray(z)
Hz = np.asarray(Hz)
sigma = np.asarray(sigma)

def chisqfunc_timesphere(arr):
	a = arr[0]
	b = arr[1]
	model = a*z+b
	chisq = np.sum( ((Hz - model)/sigma)**2 )
	return chisq
	
def E_z(z, OM, OK, OL):
    return np.sqrt((1 + z)**3 * OM + (1 + z)**2 * OK + OL)
	
def chisqfunc_lcdm(arr):
	H0 = arr[0]
	Omega_m0 = arr[1]
	Omega_l0 = arr[2]
	Omega_k0 = 0
	model = H0*E_z(z, Omega_m0, Omega_k0, Omega_l0)
	chisq = np.sum( ((Hz - model)/sigma)**2 )
	return chisq
	

Harray = np.linspace(68, 69, num=1001)
H2array = np.linspace(62, 63, num=1001)
OMarray = np.linspace(0.31, 0.33, num=1001)

tsfit = 0  # 1 if you want to fit timesphere
twoparam = 0 # 1 if you want to fit timesphere with 2 parameters
lcdmfit = 1 # 1 if you want to fit LambdaCDM


#initial_guess_timesphere = np.array([70, 70])
#initial_guess_lcdm = np.array([67.66, 0.3111, 0.6889])
#result_timesphere =  opt.minimize(chisqfunc_timesphere, initial_guess_timesphere)
#result_lcdm =  opt.minimize(chisqfunc_lcdm, initial_guess_lcdm)

if (twoparam == 1):
	res_ts = np.zeros((len(Harray), len(H2array)))
else:
	res_ts = []

res_lcdm = np.zeros((len(Harray), len(OMarray)))

for i in range(len(Harray)):

	if (tsfit == 1):

		if (twoparam == 1):
			for j in range(len(H2array)):
				res_ts[i][j] = chisqfunc_timesphere((Harray[i], H2array[j]))
		else:
			res_ts.append(chisqfunc_timesphere((Harray[i], Harray[i])))
			
	if (lcdmfit == 1):	
	
		for j in range(len(OMarray)):
			res_lcdm[i][j] = chisqfunc_lcdm((Harray[i], OMarray[j], 1-OMarray[j]))

res_ts = np.asarray(res_ts)


if (tsfit == 1):
	if (twoparam == 1):
		(inda, indb) = np.unravel_index(res_ts.argmin(), res_ts.shape)
		print(Harray[inda], H2array[indb], res_ts[inda][indb])
	else:
		print(Harray[np.argmin(res_ts)], min(res_ts))

if (lcdmfit == 1):
	(inda, indb) = np.unravel_index(res_lcdm.argmin(), res_lcdm.shape)
	print(Harray[inda], OMarray[indb], res_lcdm[inda][indb])
		


'''
fig, ax = plt.subplots(1)
plt.errorbar(z, Hz, yerr=sigma, fmt='o')
plt.plot(z, result_timesphere.x[0]*z+result_timesphere.x[1], label='Timesphere')
plt.plot(z, result_lcdm.x[0]*E_z(z, result_lcdm.x[1], 0, result_lcdm.x[2]), label='Lambda-CDM')
plt.legend(loc='upper left', prop={'size': 18})

fig.set_size_inches(8, 6, forward=True)
plt.xlabel(r'$z$', fontsize=24)
plt.ylabel(r'$H\ [\mathrm{km\ s^{-1}\ Mpc^{-1}}]$', fontsize=26)
plt.tick_params(axis='both', which='major', labelsize=18)
fig.tight_layout()

plt.savefig('Chronometer.png')
plt.show()
'''
