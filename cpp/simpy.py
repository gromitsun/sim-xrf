# simpy.py
import pyapi

import time
start = time.time()
spec = pyapi.calc(input_file = "input.txt",	output_file = "output.txt")
print time.time() - start
spec.show(xlim = [0,11], ylim = [1e-14, 1e-4])

# import matplotlib.pyplot as plt
# plt.plot(spec.y_sep.T)
# plt.plot(spec.y_vec)
# plt.yscale('log')
# plt.ylim(1e-14,1e-4)
# plt.legend(spec.labels+['total'])
# plt.show()