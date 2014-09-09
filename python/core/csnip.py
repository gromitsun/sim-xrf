import ctypes
import numpy as np
import os
libpath = os.path.dirname(os.path.realpath(__file__))
lib = ctypes.cdll.LoadLibrary(libpath+'\snip.so')

FWHM_c = lib.FWHM_c

FWHM_c.restype = ctypes.c_double # reset return types. default is c_int
FWHM_c.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double]

def FWHM(x, noise = 100., fano = 0.114):
	return FWHM_c(ctypes.c_double(x), ctypes.c_double(noise), ctypes.c_double(fano))


energy_to_channel_c = lib.energy_to_channel_c

energy_to_channel_c.restype = ctypes.c_double # reset return types. default is c_int
energy_to_channel_c.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double]

def energy_to_channel(energy, offset = 0, gain = 10):
	return energy_to_channel_c(ctypes.c_double(energy), ctypes.c_double(offset), ctypes.c_double(gain))




lsdf_c = lib.lsdf_c

lsdf_c.restype = ctypes.c_void_p # reset return types. default is c_int
lsdf_c.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.CFUNCTYPE(ctypes.c_double,ctypes.c_double), ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double]

def lsdf(E, y, FWHM = FWHM,
			f = 1.5,
			A = 75,
			M = 10,
			r = 1.3):
	
	f *= 1.
	A *= 1.
	M *= 1.
	r *= 1.
	E *= 1.
	y = 1.*y.copy()
	
	lsdf_c(E.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
			y.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
			ctypes.c_int(len(E)),
			ctypes.CFUNCTYPE(ctypes.c_double,ctypes.c_double)(FWHM),
			ctypes.c_double(f),
			ctypes.c_double(A),
			ctypes.c_double(M),
			ctypes.c_double(r))
	return y			


snip_c = lib.snip_c

snip_c.restype = ctypes.c_void_p # reset return types. default is c_int
snip_c.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.CFUNCTYPE(ctypes.c_double,ctypes.c_double),ctypes.c_double,ctypes.c_double,ctypes.c_int,ctypes.c_int,ctypes.c_double,ctypes.c_double]


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
	
	offset *= 1.
	gain *= 1.
	E *= 1.
	y = 1.*y.copy()
	
	
	snip_c(E.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
					y.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
					ctypes.c_int(len(E)),
					ctypes.CFUNCTYPE(ctypes.c_double,ctypes.c_double)(FWHM),
					ctypes.c_double(offset),
					ctypes.c_double(gain),
					ctypes.c_int(loops),
					ctypes.c_int(end_loops),
					ctypes.c_double(factor),
					ctypes.c_double(reduce_factor))
					
	return y


