#savespec
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterMathtext
import h5py
import os

# # #fonts# # #
import matplotlib
from matplotlib import rc

matplotlib.rcParams['pdf.fonttype'] = 'truetype'
fontProperties = {'family':'serif','serif':['Arial'],
    'weight' : 'normal', 'size' : '12'}
rc('font',**fontProperties)
# # #

from sim_class import *

def dim_scale(f,dset):
	ndims = len(dset.dims)
	dset.dims[ndims-2].label = 'theta'
	dset.dims[ndims-1].label = 'beta'
	dset.dims.create_scale(f['omega/theta'])
	dset.dims.create_scale(f['omega/theta_rad'])
	dset.dims.create_scale(f['omega/beta'])
	dset.dims.create_scale(f['omega/beta_rad'])
	dset.dims[ndims-2].attach_scale(f['omega/theta'])
	dset.dims[ndims-2].attach_scale(f['omega/theta_rad'])
	dset.dims[ndims-1].attach_scale(f['omega/beta'])
	dset.dims[ndims-1].attach_scale(f['omega/beta_rad'])

def savespec(rawspec, filename = 'specdata.h5', dir = os.getcwd()):
	if not os.path.exists(dir):
		os.mkdir(dir)
	rawspec_arr = flatten([rawspec])
	try:
		f = h5py.File(dir+filename, 'w')
		
		f['omega/theta'] = rawspec_arr[0].omega.theta_arr
		f['omega/theta_rad'] = rawspec_arr[0].omega.theta_rad_arr
		f['omega/beta'] = rawspec_arr[0].omega.beta_arr
		f['omega/beta_rad'] = rawspec_arr[0].omega.beta_rad_arr
		f['omega/range'] = rawspec_arr[0].omega.range
		f['omega/n_theta'] = rawspec_arr[0].omega.n_theta
		f['omega/n_beta'] = rawspec_arr[0].omega.n_beta
		
		for rawspec in rawspec_arr:
			
			group = f.require_group(rawspec.stype+'/'+rawspec.name)
			# try:
				# group = f[rawspec.stype+'/'+rawspec.name]
			# except KeyError:
				# group = f.create_group(rawspec.stype+'/'+rawspec.name)
			dset_list = rawspec.__dict__.keys()
			dset_list.remove('stype')
			dset_list.remove('name')
			dset_list.remove('omega')
			for dset in dset_list:
				try:
					group[dset] = rawspec.__getattribute__(dset)
					dim_scale(f,group[dset])
				except TypeError:
					pass
	finally:
		f.close()

def readspec(filename = 'specdata.h5', dir = os.getcwd()):
	try:
		f = h5py.File(dir+filename, 'r')
		out = []
		omega = solid_angle(angle_range=f['omega/range'].value,n_theta = f['omega/n_theta'].value,n_beta = f['omega/n_beta'].value)
		stype_list = f.keys()
		stype_list.remove('omega')
		for stype in stype_list[::-1]:
			for name in f[stype]:
				group = f[stype][name]
				out.append(spectrum(ev_mat = group['ev_mat'], intensity_mat = group['intensity_mat'], omega = omega, name = name, stype = stype))
		return out
	finally:
		f.close()

def plot_specs(specs,detector, **kwargs):
	xlim = kwargs.get('xlim',detector.channel.range)[:2]
	show_total = kwargs.get('show_total', True)
	logy = kwargs.get('logy', True)
	show = kwargs.get('show', True)
	yshift = kwargs.get('yshift', 1e-18)
	save_npz = kwargs.get('save_npz',False)
	npz_name = kwargs.get('npz_name','spec.npz')
	total = 0
	ax = plt.subplot(111)
	for x in specs:
		if (x.intensity_arr > 1e-18).any():
			plt.plot(x.ev_arr/1.e3,x.intensity_arr+yshift,label = x.name)
		else:
			plt.plot(x.ev_arr/1.e3,x.intensity_arr+yshift)
		total += x.intensity_arr
	if show_total:
		plt.plot(x.ev_arr/1.e3, total+yshift, label = 'total', linestyle='--')
		# plt.fill_between(x.ev_arr[845:880]/1.e3,1e-18,total[845:880],color = 'black',alpha = 0.15)
	ymax = total.max()
	ylim = kwargs.get('ylim',[1e-18,5*ymax])


	# Shink current axis by 20%
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
	# plt.legend(loc=0,ncol=3)
	plt.legend(ncol=1,bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	if logy:
		plt.yscale('log')
	plt.xlim(xlim)
	plt.ylim(ylim)
	plt.xlabel('E (KeV)')
	plt.ylabel(r'$I(E)/I_0(E_0)$')
	if show:
		plt.show()
	if save_npz:
		np.savez(npz_name,ev_arr = x.ev_arr,total = total)
		