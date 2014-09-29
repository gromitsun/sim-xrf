# pyapi.py
# # # This file must be put in the same directory as main.cpp.
from python.elementlookup import number2symbol # Must be imported before loading the DLL.
import ctypes
import numpy as np
import os

libpath = os.path.dirname(os.path.realpath(__file__))+'/Lib'
try:
	lib = ctypes.cdll.LoadLibrary(libpath+'/libsim.so')
except:
	lib = ctypes.cdll.LoadLibrary(libpath+'/libsim.dll')

from python.spectrum import *


lib.sim.restype = None
lib.sim.argtypes = [ctypes.c_char_p, 
				ctypes.c_char_p, 
				ctypes.POINTER(ctypes.c_double), 
				ctypes.POINTER(ctypes.c_double), 
				ctypes.POINTER(ctypes.c_int), 
				ctypes.POINTER(ctypes.c_int), 
				ctypes.POINTER(ctypes.c_int), 
				ctypes.POINTER(ctypes.c_double), 
				ctypes.POINTER(ctypes.c_double), 
				ctypes.POINTER(ctypes.c_double), 
				ctypes.POINTER(ctypes.c_double), 
				ctypes.POINTER(ctypes.c_double), 
				ctypes.POINTER(ctypes.c_double), 
				ctypes.POINTER(ctypes.c_int), 
				ctypes.POINTER(ctypes.c_int), 
				ctypes.POINTER(ctypes.c_double), 
				ctypes.POINTER(ctypes.c_double), 
				ctypes.POINTER(ctypes.c_int)]

def sim(input_file, 
		output_file, 
		y_vec,
		y_sep,
		Z_vec,
		row,
		lines,
		xrf_ev,
		xrf_y,
		comp_ev,
		comp_y,
		ray_y,
		det,
		n_channels, 
		win_mat, 
		il,
		sa,
		nout):
		
	if not os.path.isfile(input_file):
		raise IOError("File %s does not exit!" % input_file)
	
	lib.sim(ctypes.c_char_p(input_file), 
		ctypes.c_char_p(output_file), 
		y_vec.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
		y_sep.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
		Z_vec.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
		row.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
		lines.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
		xrf_ev.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
		xrf_y.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
		comp_ev.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
		comp_y.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
		ray_y.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
		det.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
		n_channels.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
		win_mat.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
		il.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
		sa.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
		nout.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))


def calc(input_file = "input.txt",	
	output_file = "output.txt", 
	nout = [1000, 30, 500, 500]): # N of channels, N of Z, N of lines, N of thetas
	y_vec = np.zeros(nout[0])
	y_sep = np.zeros(nout[0]*(nout[1]+2))
	Z_vec = np.zeros(nout[1], dtype = int)
	row = np.zeros(nout[1]+2, dtype = int)
	lines = np.zeros(nout[2], dtype = int)
	xrf_ev = np.zeros(nout[2])
	xrf_y = np.zeros(nout[2])
	comp_ev = np.zeros(nout[3])
	comp_y = np.zeros(nout[3])
	ray_y = np.zeros(1)
	det = np.zeros(9)
	n_channels = np.zeros(1, dtype = int)
	win_mat = np.zeros(20, dtype = int)
	il = np.zeros(3)
	sa = np.zeros(6)
	nout = np.array(nout, dtype = int)
		
	sim(input_file, 
			output_file,
			y_vec,
			y_sep,
			Z_vec,
			row,
			lines,
			xrf_ev,
			xrf_y,
			comp_ev,
			comp_y,
			ray_y,
			det,
			n_channels,
			win_mat,
			il,
			sa,
			nout)
	
	xrf = Xrf(None, xrf_y[:nout[2]], xrf_ev[:nout[2]], lines[:nout[2]], Z_vec[:nout[1]], row[:nout[1]+1])
	comp = Compton(None, comp_y[:nout[3]], comp_ev[:nout[3]])
	ray = Rayleigh(None, ray_y, il[0])
	win_mat = ''.join(chr(ch) for ch in win_mat)
	_det = Detector(Channel(*det[:2], n_channels = n_channels),
					Response(*det[2:-2], ev_gain = det[1]),
					Window(win_mat, *det[-2:]))
	_il = Illumination(*il)
	omega = solid_angle(sa[:4], *sa[-2:])
	
	# labels = Z_vec[:nout[1]].tolist()
	labels = [number2symbol(Z) for Z in Z_vec[:nout[1]]]
	labels.append('Rayleigh')
	labels.append('Compton')
	
	return Spectrum(y_vec[:n_channels], y_sep[:n_channels*(nout[1]+2)].reshape(-1, n_channels), labels, xrf, ray, comp, _det, _il, omega)
	

	# print y_vec
	# print Z_vec
	# print row
	# print lines
	# print xrf_ev
	# print xrf_y
	# print comp_ev
	# print comp_y
	# print ray_y
	# print det
	# print n_channels
	# print win_mat
	# print il
	# print sa
	# print nout
	
# spec = calc()

# print spec.detector.window.material
# print spec.labels
# print spec.y_sep[0]