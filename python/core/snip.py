import numpy as np
from scipy.optimize import curve_fit

def FWHM(x, noise = 100, fano = 0.114):
	sigma = np.sqrt((noise/2.3548)**2+3.58*fano*x)
	return 2.3548*sigma

def fit_FWHM(x, F):
	def _FWHM(x, noise, fano):
		return (noise/2.3548)**2+3.58*fano*x
	popt, pcov = curve_fit(_FWHM, x, (F/2.3548)**2, p0 = [100,0.114])
	return popt


def energy_to_channel(energy, offset = 2.97, gain = 12.26952):
	return 1.*(energy - offset)/gain



# # # Low statistics digital filter

def lsdf(E, y, FWHM = FWHM,
			f = 1.5,
			A = 75,
			M = 10,
			r = 1.3):
	def reduce(x, length_start):
		for i in range(length_start):
			length = length_start - i
			if x < length:
				raise IndexError
			L = y[x-length:x].sum()
			R = y[x+1:x+length+1].sum()
			S = y[x] + L + R
			slope = (R + 1.)/(L + 1.)
			if S < M or S < A*np.sqrt(y[x]) or (slope >= 1./r and slope <= r):
				return S/(2.*length+1)
		print 'Not found for x = %d!' % x
		return y[x]
	
	y_out = y.copy()
	
	for x in range(len(E)):
		try:
			len_0 = int(energy_to_channel(f*FWHM(E[x]), E[0], E[1]-E[0]))
			y_out[x] = reduce(x, len_0)
		except IndexError:
			pass
	return y_out

# # # Peak-clipping

def snip(E,y,FWHM = FWHM, offset = 0., gain = 10., **kwargs):
	det = kwargs.get('detector')
	loops = kwargs.get('loops',24)
	end_loops = kwargs.get('end_loops',8)
	reduce_factor = kwargs.get('reduce_factor',np.sqrt(2))
	factor = kwargs.get('factor',2)
	
	
	if det != None:
		FWHM = det.response.FWHM
		offset = det.channel.offset
		gain = det.channel.gain
		
	def G(y):
		return np.log(np.log(y+1)+1)
		
	def w(x, factor = 2):
		return energy_to_channel(factor*FWHM(E[x]), offset = offset, gain = gain)

	def G_inv(z):
		return np.exp(np.exp(z)-1)-1
	
	z_out = G(y)
	for i in range(loops):
		if i >= loops - end_loops:
			factor /= 1.*reduce_factor
		z = z_out.copy()
		for x in range(len(E)):
			try:
				_w = w(x, factor = factor)
				if _w > x:
					raise IndexError
				z_bar = (z[x+_w]+z[x-_w])/2.
				z_out[x] = min(z[x],z_bar)
			except IndexError:
				pass
	return G_inv(z_out)


