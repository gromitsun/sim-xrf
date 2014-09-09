import ctypes
import numpy as np
import os
libpath = os.path.dirname(os.path.realpath(__file__))
lib = ctypes.cdll.LoadLibrary(libpath+'\libscatt_bg.so')

scatt_bg_c = lib.scatt_bg

scatt_bg_c.restype = ctypes.c_void_p # reset return types. default is c_int
scatt_bg_c.argtypes = [ctypes.c_double, ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_int]

subtend_c = lib.subtend

subtend_c.restype = ctypes.c_double # reset return types. default is c_int
subtend_c.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double]


def scatt_bg(kev, Z_max = 98, theta_max = 90):
	# kev = np.asarray(kev)
	out = np.zeros(Z_max*theta_max)
	# Z_max = np.asarray(Z_max)
	# theta_max = np.asarray(theta_max)
	scatt_bg_c(ctypes.c_double(kev),out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),ctypes.c_int(Z_max),ctypes.c_int(theta_max))
	return out
	
def subtend(theta0,theta1,beta0,beta1):
	return subtend_c(c_double(np.radians(theta0)),c_double(np.radians(theta1)),c_double(np.radians(beta0)),c_double(np.radians(beta1)))
	
