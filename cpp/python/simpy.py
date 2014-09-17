import ctypes
import numpy as np
import os
libpath = os.path.dirname(os.path.realpath(__file__))+'/../Lib'
try:
	lib = ctypes.cdll.LoadLibrary(libpath+'/libsim.so')
except WindowsError:
	lib = ctypes.cdll.LoadLibrary(libpath+'/libsim.dll')
	
lib.sim.restype = ctypes.c_void_p
lib.sim.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]

def sim(input_file, 
		output_file, 
		y_vec,
		Z_vec,
		row,
		lines,
		xrf_ev,
		xrf_y,
		comp_ev,
		comp_y,
		ray_y,
		ev_offset,
		ev_gain,
		n_channels, 
		nout):
	lib.sim(ctypes.c_char_p(input_file), 
		ctypes.c_char_p(output_file), 
		y_vec.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
		Z_vec.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
		row.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
		lines.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
		xrf_ev.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
		xrf_y.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
		comp_ev.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
		comp_y.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
		ray_y.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
		ev_offset.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
		ev_gain.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
		n_channels.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
		nout.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))


def get_data(input_file = "input.txt",	output_file = "output.txt"):
	pass


input_file = "input.txt" 
output_file = "output.txt"
y_vec = np.zeros(2100)
Z_vec = np.zeros(100, dtype = int)
row = np.zeros(100, dtype = int)
lines = np.zeros(100, dtype = int)
xrf_ev = np.zeros(100)
xrf_y = np.zeros(100)
comp_ev = np.zeros(30)
comp_y = np.zeros(30)
ray_y = np.zeros(1)
ev_offset = np.zeros(1)
ev_gain = np.zeros(1)
n_channels = np.zeros(1, dtype = int)
nout = np.zeros(8, dtype = int)
		
sim(input_file, 
		output_file,
		y_vec,
		Z_vec,
		row,
		lines,
		xrf_ev,
		xrf_y,
		comp_ev,
		comp_y,
		ray_y,
		ev_offset,
		ev_gain,
		n_channels,
		nout)

print y_vec
print Z_vec
print row
print lines
print xrf_ev
print xrf_y
print comp_ev
print comp_y
print ray_y
print ev_offset
print ev_gain
print n_channels
print nout