# genspec
from sim_common import *

def _mat_reduce(ev_mat, intensity_mat):
	sum_axes =  tuple(np.where(np.array(ev_mat.shape) - np.array(intensity_mat.shape))[0])
	_intensity_mat = np.sum(intensity_mat,axis=sum_axes,keepdims=True)
	return ev_mat, _intensity_mat

def genspec_raw(ev_mat, intensity_mat):
	"""
	generate spectrum from intensity_mat, w/o detector channel bins
	"""
	ev_arr, intensity_arr = (x.flatten() for x in _mat_reduce(ev_mat, intensity_mat))
	return ev_arr, intensity_arr

def genspec_bin(ev_arr, intensity_arr, energy_range,out='center'):
	"""
	generate spectrum without detector response
	out: output format for energy.
		center:	return the centers of the energy bins
		bins:	return the energy bin edges, including the rightmost
	"""
	energy_range_1 = np.array(energy_range)+ [0, energy_range[2], 0]
	energy_bins = np.arange(*energy_range_1) - energy_range[-1]/2. #energy bin edges, including the rightmost
	spec, energy = np.histogram(ev_arr, weights = intensity_arr, bins = energy_bins)
	if out == 'center':
		energy = (energy[:-1]+energy[1:])/2.
	return energy, spec

def genspec_det(ev_arr, intensity_arr, detector):
	"""
	generate spectrum with detector response
	"""
	# process inputs
	ev_arr = np.array(ev_arr).flatten()
	intensity_arr = np.array(intensity_arr).flatten()
	# Gaussian function
	noise = detector.response.noise
	fano = detector.response.fano
	gamma = detector.response.gamma
	fs = detector.response.fs
	ft = detector.response.ft
	gain = detector.response.gain
	energy_range = detector.channel.range
	ev_ch_arr = detector.channel.centers
	
	G = np.zeros(np.shape(ev_ch_arr))
	S = np.zeros(np.shape(ev_ch_arr))
	T = np.zeros(np.shape(ev_ch_arr))

	# reshape spectrum arr
	ev_arr = np.array(ev_arr.flatten(),ndmin=2)
	intensity_arr = np.array(intensity_arr.flatten(),ndmin=2)
	ev_ch_arr = np.array(ev_ch_arr,ndmin=2).T
	
	# apply window filter
	if detector.window != None:
		intensity_arr = intensity_arr*detector.window.transmission(ev_arr)
	
	# apply detector response function
	sigma = np.sqrt((noise/2.3548)**2+3.58*fano*ev_arr)
	G = np.sum(intensity_arr*gain/(sigma*np.sqrt(2*np.pi))*np.exp(-(ev_arr-ev_ch_arr)**2/(2*sigma**2)),axis=1)
	S = np.sum(intensity_arr*gain/(2.*ev_arr)*erfc((ev_ch_arr-ev_arr)/(np.sqrt(2)*sigma)),axis=1)
	T = np.sum(intensity_arr*gain/(2.*gamma*sigma*np.exp(-1./(2*gamma**2)))*np.exp((ev_ch_arr-ev_arr)/(gamma*sigma))*erfc((ev_ch_arr-ev_arr)/(np.sqrt(2)*sigma)+1/(np.sqrt(2)*gamma)),axis=1)
	
	return ev_ch_arr.flatten(), G+fs*S+ft*T

