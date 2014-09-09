import sys, os
sys.path.append(os.path.dirname(os.path.realpath(__file__))+'/../core/')
from sim_class import *



angle_range_gen1 = lambda x: [90-x, 90+x, -x, x]

angle_range_gen = lambda x: [0, x, 0, 360]


y = np.array([solid_angle(angle_range = angle_range_gen(x)).subtend for x in range(0,91)])
y1 = np.array([solid_angle(angle_range = angle_range_gen1(x)).subtend for x in range(0,91)])

import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rc

matplotlib.rcParams['pdf.fonttype'] = 'truetype'
fontProperties = {'family':'serif','serif':['Arial'],
    'weight' : 'normal', 'size' : '12'}
rc('font',**fontProperties)

plt.plot(y,label = r'$0^\circ$ and $180^\circ$')
plt.plot(y1,'r',label = r'$90^\circ$')
plt.xlabel(r'Collection semi-angle $\phi$ (deg)')
plt.ylabel('Total subtended solid angle (sr)')
plt.legend(loc = 'upper left')
plt.grid('on')
plt.show()