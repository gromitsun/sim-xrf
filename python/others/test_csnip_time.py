from csnip import *
from time import time

s = time()
a = np.arange(10e3)
b = a*1.2
snip(a,b,FWHM = lambda x: 1,loops = 24,end_loops = 8,gain=1.)

print time()-s