# spectrum.py
from geometry import *
from detector import *

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
		total = None,
		y_vec = None,
		labels = [],
		xrf = Xrf(),
		ray = Rayleigh(),
		comp = Compton(),
		det = Detector(),
		il = Illumination(),
		omega = solid_angle()):
		self.total = total
		self.y_vec = y_vec
		self.labels = labels
		self.xrf = xrf
		self.rayleigh = ray
		self.compton = comp
		self.detector = det
		self.illumination = il
		self.omega = omega
	