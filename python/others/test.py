import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scatt_bg_c import *
im = scatt_bg(10,5,5).reshape(5,5)

np.savez('scatt_10kev.npz',im=im)
im=np.load('scatt_10kev.npz')['im']
print im
plt.imshow(im,extent=[1,90,1,98],norm = LogNorm())
plt.colorbar()
plt.show()
