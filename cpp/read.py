import numpy as np
import matplotlib.pyplot as plt
a = np.genfromtxt('out.txt',invalid_raise = False,skiprows=3)
plt.plot(a)
plt.yscale('log')
plt.ylim(1e-10,1e-5)
plt.xlim(0,1300)

plt.show()