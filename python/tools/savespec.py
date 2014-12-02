#savespec
import h5py
import os

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