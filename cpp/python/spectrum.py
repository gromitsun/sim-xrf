# spectrum.py
from geometry import *
from detector import *

import matplotlib.pyplot as plt
# from matplotlib.colors import LogNorm
# from matplotlib.ticker import LogFormatterMathtext

# # # #fonts# # #
# import matplotlib
# from matplotlib import rc

# matplotlib.rcParams['pdf.fonttype'] = 'truetype'
# fontProperties = {'family':'serif','serif':['Arial'],
    # 'weight' : 'normal', 'size' : '12'}
# rc('font',**fontProperties)
# # # #


class Xrf(object):
	def __init__(self, y_mat = None,
				y_vec = None, ev_vec = None, lines = None, Z_vec = None, row = None):
		self.y_mat = y_mat
		self.y_vec = y_vec
		self.ev_vec = ev_vec
		self.lines = lines
		self.Z_vec = Z_vec
		self.row = row
		
class Rayleigh(object):
	def __init__(self, y_mat = None, 
				y = None, ev0 = None):
		self.y_mat = y_mat
		self.y = y
		self.ev0 = ev0

		
class Compton(object):
	def __init__(self, y_mat = None, 
				y_vec = None, ev_vec = None):
		self.y_mat = y_mat
		self.y_vec = y_vec
		self.ev_vec = ev_vec
		

class Spectrum():
	def __init__(self,
		y_vec = None,
		y_sep = None,
		labels = [],
		xrf = Xrf(),
		ray = Rayleigh(),
		comp = Compton(),
		det = Detector(),
		il = Illumination(),
		omega = solid_angle()):
		self.y_vec = y_vec
		self.y_sep = y_sep
		self.labels = labels
		self.xrf = xrf
		self.rayleigh = ray
		self.compton = comp
		self.detector = det
		self.illumination = il
		self.omega = omega
	def show(self, **kwargs):
		xlim = kwargs.get('xlim')
		show_total = kwargs.get('show_total', True)
		logy = kwargs.get('logy', True)
		show = kwargs.get('show', True)
		yshift = kwargs.get('yshift', 0)
		ax = plt.subplot(111)
		plt.plot(self.detector.channel.ev_arr/1.e3, self.y_sep.T+yshift)
		labels = self.labels
		if show_total:
			plt.plot(self.detector.channel.ev_arr/1.e3, self.y_vec+yshift, linestyle='--')
			# plt.fill_between(x.ev_arr[845:880]/1.e3,1e-18,total[845:880],color = 'black',alpha = 0.15)
			labels += ['Total']
		ymax = self.y_vec.max()
		ylim = kwargs.get('ylim',[1e-18,5*ymax])


		# Shink current axis by 20%
		box = ax.get_position()
		ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
		# plt.legend(loc=0,ncol=3)
		plt.legend(self.labels, ncol=1,bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
		if logy:
			plt.yscale('log')
		plt.xlim(xlim)
		plt.ylim(ylim)
		plt.xlabel('E (KeV)')
		plt.ylabel(r'$I(E)/I_0(E_0)$')
		if show:
			plt.show()
	
	def save_npz(self, filename = 'spec.npz'):
		np.savez(filename, ev_arr = self.detector.channel.ev_arr, y_vec = self.y_vec, y_sep = self.y_sep)
		
	def save_hdf(self, filename = 'spec.hdf5'):
		pass
