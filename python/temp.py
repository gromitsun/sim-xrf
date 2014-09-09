from sim_class import *



angle_range_gen1 = lambda x: [90-x, 90+x, -x, x]

angle_range_gen = lambda x: [0, x, 0, 360]


y = np.array([solid_angle(angle_range = angle_range_gen(x)).subtend for x in range(0,91)])
y1 = np.array([solid_angle(angle_range = angle_range_gen1(x)).subtend for x in range(0,91)])

import matplotlib.pyplot as plt
plt.plot(y,label = r'$0^\circ$ and $180^\circ$')
plt.plot(y1,'r',label = r'$90^\circ$')
plt.xlabel('Opening angle (deg)')
plt.ylabel('Total solid angle (sr)')
plt.legend(loc = 'upper left')
plt.grid('on')
plt.show()