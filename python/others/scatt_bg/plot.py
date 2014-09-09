import numpy as np
import matplotlib.pyplot as plt
data = np.load("30kev.npz")
ray = data['ray'][:98]
comp = data['comp'][:98]
im = ray+comp
im = im[:,6:] - im[:,5:-1]
neg = im.copy()
pos = im.copy()
neg[np.where(neg>=0)] = np.nan
pos[np.where(pos<=0)] = np.nan
plt.imshow(im,origin='lower',cmap='bwr',vmin=im.min(),vmax=-im.min())
# plt.imshow(neg*-1,origin='lower',cmap='Blues')
# plt.colorbar()
# plt.imshow(pos,origin='lower',cmap='Reds')
plt.colorbar()
plt.show()