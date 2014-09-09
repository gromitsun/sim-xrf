# plot snip evolutions
import sys, os
sys.path.append(os.path.realpath('./core/'))
import numpy as np
import matplotlib.pyplot as plt
import snip

npz_file = 'spec_100um_ca50.npz'

E = np.load(npz_file)['ev_arr']
y = np.load(npz_file)['total']
y1 = y.copy()
n1,n2 = np.where(E==8130)[0][0],np.where(E==9130)[0][0]

E = E[n1:n2]
y = y[n1:n2]
y1 = y1[n1:n2]


for x in range(16):
	y = snip.snip(E,y,loops = 1,end_loops = 0)
	plt.figure(str(x+1)+'_0')
	plt.title('standard loops: '+str(x+1)+', ending loops: 0')
	plt.plot(E,y,'g--',label='SNIP')
	plt.plot(E,y1,'r',label='spectrum')
	# plt.yscale('log')
	plt.ylim(2e-8,2.5e-8)
	plt.xlim(8400,8800)
	plt.savefig(str(x+1)+'_0.pdf',format = 'pdf')
	
factor = 2
for x in range(8):
	factor /= np.sqrt(2)
	y = snip.snip(E,y,loops = 1,end_loops = 0, factor = factor)
	plt.figure('16_'+str(x+1))
	plt.title('standard loops: 0, ending loops: '+str(x+1))
	plt.plot(E,y,'g--',label='SNIP')
	plt.plot(E,y1,'r',label='spectrum')
	# plt.yscale('log')
	plt.ylim(2e-8,2.5e-8)
	plt.xlim(8400,8800)
	plt.savefig('16_'+str(x+1)+'.pdf',format = 'pdf')
# plt.show()