# detector.py
import numpy as np
import xraylib as xrl

def energy_to_channel(energy, offset = 0., gain = 10.):
	return (energy - offset)//gain
	
class Channel(object):
	def __init__(self, ev_offset = 0, ev_gain = 10, n_channels = 2048):
		self.ev_offset = ev_offset
		self.ev_gain = ev_gain
		self.n_channels = n_channels
		self.ev_arr = np.arange(ev_offset, ev_offset+ev_gain*(n_channels+1), ev_gain)
	def channel(self, ev):
		return energy_to_channel(ev, self.ev_offset, self.ev_gain)

class Response(object):
	def __init__(self,
		noise = 100,
		fano = 0.114,
		gamma = 2.5,
		fs = 0.03,
		ft = 0.02,
		ev_gain = None):
		self.noise = noise
		self.fano = fano
		self.gamma = gamma
		self.fs = fs
		self.ft = ft
		self.ev_gain = ev_gain
	def FWHM(self, x, **kwargs):
		noise = kwargs.get('noise',self.noise)
		fano = kwargs.get('fano',self.fano)
		sigma = np.sqrt((noise/2.3548)**2+3.58*fano*x)
		return 2.3548*sigma

		
class Window(object):
	def __init__(self,
		material = 'Be',
		thickness = 24e-4,
		density = None):
		self.material = material
		self.thickness = thickness
		if density == None:
			density = xrl.ElementDensity(xrl.SymbolToAtomicNumber(material))
			if not density:
				density = 1
		self.density = density
	def transmission(self, ev):
		_mac = xrl.CS_Total_CP(material, ev/1000.)
		# _mac = compound(CP=self.material).mac_total(ev)
		return np.exp(-_mac*self.density*self.thickness)
		
class Detector(object):
	def __init__(self, channel = Channel(), response = Response(), window = Window()):
		self.channel = channel
		self.response = response
		response.ev_gain = channel.ev_gain
		self.window = window