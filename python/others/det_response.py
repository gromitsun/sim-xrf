# plot of detector response function
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf
def erfc(x):
	return 1 - erf(x)
	
Ej = 8e3

GAIN = 10
NOISE = 100
FANO = 0.114
gamma = 2.5
fs = 0.03
ft = 0.02
def sigma(Ej):
	return np.sqrt((NOISE/2.3548)**2+3.58*FANO*Ej)

def G(Ei,Ej):
	return GAIN/(sigma(Ej)*np.sqrt(2*np.pi))*np.exp(-(Ej-Ei)**2/(2*sigma(Ej)**2))
	
def S(Ei,Ej):
	return GAIN/(sigma(Ej)*np.sqrt(2*Ej))*erfc((Ei-Ej)/(sigma(Ej)*np.sqrt(2)))
	
def T(Ei,Ej):
	return GAIN/(2*gamma*sigma(Ej)*np.exp(-1/(2*gamma**2)))*np.exp((Ei-Ej)/(gamma*sigma(Ej)))*erfc((Ei-Ej)/(sigma(Ej)*np.sqrt(2))+1/(gamma*np.sqrt(2)))
	
def P(Ei,Ej):
	return G(Ei,Ej)+fs*S(Ei,Ej)+ft*T(Ei,Ej)
	
Ei = np.arange(0,10e3,10)


import matplotlib
from matplotlib import rc

matplotlib.rcParams['pdf.fonttype'] = 'truetype'
fontProperties = {'family':'serif','serif':['Arial'],
    'weight' : 'normal', 'size' : '12'}
rc('font',**fontProperties)



plt.plot(Ei/1.e3, G(Ei,Ej)+1e-30,'b-',label = 'G(Ei,Ej)')
plt.plot(Ei/1.e3, fs*S(Ei,Ej)+1e-30,'y-',label = 'fs*S(Ei,Ej)')
plt.plot(Ei/1.e3, ft*T(Ei,Ej)+1e-30,'r-',label = 'ft*T(Ei,Ej)')
plt.plot(Ei/1.e3, P(Ei,Ej)+1e-30,'g-',label = 'P(Ei,Ej)')
plt.legend(loc='lower left')
plt.xlabel('Energy (KeV)')
plt.ylabel('Intensity')
plt.yscale('log')
plt.yticks(np.logspace(-14,0,8))
plt.ylim(1e-15,1)


plt.show()