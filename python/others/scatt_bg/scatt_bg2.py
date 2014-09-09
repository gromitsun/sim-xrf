import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys, os
sys.path.append(os.path.dirname(os.path.realpath(__file__))+'/../../core/')
from sim_calcspec import *

protein = compound(CP = 'C30H50O10N9S')
water = compound(CP = 'H2O')
Zn = compound(CP = 'Zn')

p1 = mix_compound([water, protein],[.75,.15])
p2 = compound(CP = 'Si')
p1=p2

ev = 1.2e4
theta = np.radians(np.arange(1,181))
beta = np.radians(np.arange(-90,91))

theta, beta = np.meshgrid(theta,beta)

scatt = p1.dmac_rayleigh_pol(ev,theta,beta)+p1.dmac_compton_pol(ev,theta,beta)


plt.imshow(scatt,origin = 'lower',extent=[1,180,-90,90],norm = LogNorm())
plt.colorbar()
plt.xlabel(r'Scattering angle $\theta$ (deg)') 
plt.ylabel(r'Polarization angle $\beta$ (deg)') 
plt.show()