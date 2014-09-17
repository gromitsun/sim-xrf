import numpy as np
import matplotlib.pyplot as plt
a = np.genfromtxt('out.txt',invalid_raise = False,skiprows=0)
plt.plot(a.T)
# plt.plot(np.sum(a[1:], axis = 0), 'ro')
plt.yscale('log')
plt.ylim(1e-20,1e-4)
plt.xlim(0,1300)

plt.show()