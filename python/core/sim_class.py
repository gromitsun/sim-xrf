from sim_cs_func import *

class solid_angle(object):
	def __init__(self, angle_range = [0, 180, -90, 90], **kwargs):
		self.range = angle_range
		self.n_theta = kwargs.get('n_theta', (angle_range[1]-angle_range[0])*2+1)
		self.n_beta = kwargs.get('n_beta', (angle_range[3]-angle_range[2])*2+1)
		self.theta_arr = np.linspace(angle_range[0], angle_range[1], self.n_theta)
		self.beta_arr = np.linspace(angle_range[2], angle_range[3], self.n_beta)
		self.theta_inc = self.theta_arr[min(1,self.n_theta-1)] - self.theta_arr[0]
		self.beta_inc = self.beta_arr[min(1,self.n_beta-1)] - self.beta_arr[0]
		self.theta_rad_arr = deg2rad(self.theta_arr)
		self.beta_rad_arr = deg2rad(self.beta_arr)
		self.theta_mesh, self.beta_mesh = np.meshgrid(self.theta_rad_arr, self.beta_rad_arr)
		self.delta_theta = deg2rad(np.abs(self.theta_inc))
		self.delta_beta = deg2rad(np.abs(self.beta_inc))
		self.center = [(angle_range[1] + angle_range[0])/2.,(angle_range[2] + angle_range[3])/2.]
		self.subtend_calc = np.sum(np.sin(self.theta_mesh)*self.delta_theta*self.delta_beta)
		self.subtend = (np.cos(self.theta_rad_arr[0])-np.cos(self.theta_rad_arr[-1]))*(self.beta_rad_arr[-1]-self.beta_rad_arr[0])
	def shift(self, new_center):
		pass
	def slice(self,angle_range=[0,180,-180,180]):
		theta_slice = np.logical_and(self.theta_arr>=angle_range[0],self.theta_arr<=angle_range[1])
		beta_slice = np.logical_and(self.beta_arr>=angle_range[2],self.beta_arr<=angle_range[3])
		omega_slice = np.multiply(*np.meshgrid(theta_slice,beta_slice))
		return omega_slice
	def domega(self,theta_rad):
		return (np.cos(theta_rad-self.delta_theta/2.)-np.cos(theta_rad+self.delta_theta/2.))*self.delta_beta

		
class compound(object):
	def __init__(self,Z_arr=None,p_arr=None,CP=None):
		if CP != None:
			Z_arr, p_arr = compound_parser(CP)
		self.Z_arr = symbol2number(Z_arr)
		self.Z_arr_symbol = number2symbol(Z_arr)
		self.p_arr = p_arr
		self.molecular_weight = molecular_weight(Z_arr,p_arr)
	def _mac_total(self,ev):
		return np.sum(mac_total(ev,self.Z_arr)*self.p_arr*atomic_weight(self.Z_arr))/self.molecular_weight
	def _dmac_rayleigh_pol(self,ev,theta_rad,beta_rad):
		return np.sum(dcs_rayleigh_pol(ev,self.Z_arr,theta_rad,beta_rad)*self.p_arr*N_A/self.molecular_weight)
	def _dmac_compton_pol(self,ev,theta_rad,beta_rad):
		return np.sum(dcs_compton_pol(ev,self.Z_arr,theta_rad,beta_rad)*self.p_arr*N_A/self.molecular_weight)
	def mac_total(self,ev):
		return np.array(np.frompyfunc(self._mac_total,1,1)(ev),dtype=np.float64)
	def dmac_rayleigh_pol(self,ev,theta_rad,beta_rad):
		return np.array(np.frompyfunc(self._dmac_rayleigh_pol,3,1)(ev,theta_rad,beta_rad),dtype=np.float64)
	def dmac_compton_pol(self,ev,theta_rad,beta_rad):
		return np.array(np.frompyfunc(self._dmac_compton_pol,3,1)(ev,theta_rad,beta_rad),dtype=np.float64)
	
def mix_compound(compound_arr,weight_ratio_arr):
	Z_arr = np.array([],dtype = np.int32)
	p_arr = np.array([])
	for i in range(len(compound_arr)):
		Z_arr_1 = compound_arr[i].Z_arr
		p_arr_1 = np.array(compound_arr[i].p_arr)*weight_ratio_arr[i]/compound_arr[i].molecular_weight
		for j in range(len(Z_arr_1)):
			if Z_arr_1[j] in Z_arr:
				p_arr[np.where(Z_arr == Z_arr_1[j])] += p_arr_1[j]
			else:
				Z_arr = np.append(Z_arr, Z_arr_1[j])
				p_arr = np.append(p_arr, p_arr_1[j])
	return compound(Z_arr,p_arr)
	
class sample(object):
	def __init__(self,compound,density,thickness,layer=None):
		self.compound = compound
		self.density = density
		self.thickness = thickness
		self.layer = layer

class sample_multilayer(object):
	def __init__(self,sample_arr):
		self.nlayers = len(sample_arr)
		self.layers = np.arange(self.nlayers)
		self.sample_arr = sample_arr
		self.density_arr = np.array([x.density for x in sample_arr])
		self.thickness_arr = np.array([x.thickness for x in sample_arr])
		for i in range(len(sample_arr)):
			sample_arr[i].layer = i
		self.thickness_total = np.sum(self.thickness_arr)

		
		
		
class illumination(object):
	def __init__(self, ev0, psi, alpha=0):
		self.ev0 = ev0
		self.kev0 = ev0/1.e3
		self.psi = psi
		self.alpha = alpha
		self.psi_rad = deg2rad(self.psi)
		self.alpha_rad = deg2rad(self.alpha)

def energy_to_channel(energy, offset = 0., gain = 10.):
	return 1.*(energy - offset)/gain		

class channel(object):
	def __init__(self, ev_offset = 0, ev_gain = 10, n_channels = 2048):
		self.offset = ev_offset
		self.gain = ev_gain
		self.n_channels = n_channels
 		self.centers = np.arange(ev_offset, ev_offset+ev_gain*n_channels, ev_gain)
		self.ev_high = self.centers[-1]
		self.range = [ev_offset, self.ev_high, ev_gain]
		self.bins = np.arange(ev_offset, ev_offset+ev_gain*(n_channels+1), ev_gain)-ev_gain/2.
	def channel(self, ev):
		return energy_to_channel(ev, self.offset, self.gain)

class response(object):
	def __init__(self,
		noise = 100,
		fano = 0.114,
		gamma = 2.5,
		fs = 0.03,
		ft = 0.02,
		gain = None):
		self.noise = noise
		self.fano = fano
		self.gamma = gamma
		self.fs = fs
		self.ft = ft
		self.gain = gain
	def FWHM(self, x, **kwargs):
		noise = kwargs.get('noise',self.noise)
		fano = kwargs.get('fano',self.fano)
		sigma = np.sqrt((noise/2.3548)**2+3.58*fano*x)
		return 2.3548*sigma

		
class window(object):
	def __init__(self,
		material = 'Be',
		thickness = 24e-4,
		density = None):
		self.material = material
		self.thickness = thickness
		if density == None:
			density = element_density(symbol2number(material))
			if not density:
				density = 1
		self.density = density
	def transmission(self, ev):
		# _mac = CS_Total_CP(material, ev/1000.)
		_mac = compound(CP=self.material).mac_total(ev)
		return np.exp(-_mac*self.density*self.thickness)
		
class detector(object):
	def __init__(self, channel = channel(), response = response(), omega=solid_angle(), window = window()):
		self.channel = channel
		self.response = response
		response.gain = channel.gain
		self.omega = omega
		self.window = window

		
class spectrum(object):
	def __init__(self, **kwargs):
		self.ev_arr = np.array(kwargs.get('ev_arr'))
		self.intensity_arr = np.array(kwargs.get('intensity_arr'))
		self.ev_mat = np.array(kwargs.get('ev_mat'))
		self.intensity_mat = np.array(kwargs.get('intensity_mat'))
		self.omega = kwargs.get('omega',solid_angle())
		self.name = kwargs.get('name')
		self.stype = kwargs.get('stype')

def sum_spec(rawspec_arr,*args,**kwargs):
	mat = kwargs.get('mat')
	arr = kwargs.get('arr')
	# read in args
	if args:
		rawspec_arr = flatten([rawspec_arr]+list(args))
	name = kwargs.get('name',rawspec_arr[0].name)
	stype = kwargs.get('stype',rawspec_arr[0].stype)
	out = spectrum(name=name,stype=stype)
	if mat:
		_intensity_mat = np.concatenate([[rawspec.intensity_mat] for rawspec in rawspec_arr],axis=0)
		_intensity_mat = np.sum(_intensity_mat, axis = 0)
		out.ev_mat = rawspec.ev_mat
		out.intensity_mat = _intensity_mat
	if arr:
		_intensity_arr = np.array([rawspec.intensity_arr for rawspec in rawspec_arr])
		_intensity_arr = np.sum(_intensity_arr, axis = 0)
		out.ev_arr = rawspec.ev_arr
		out.intensity_arr = _intensity_arr
	return out