from sim_class import *
from scipy.integrate import dblquad

def subtend(angle_range):
	return (np.cos(angle_range[0])-np.cos(angle_range[1]))*(angle_range[3]-angle_range[2])

def int_spec(diff_spec, angle_range, *args, **kwargs):
	ar_rad = np.radians(angle_range)
	a, b = ar_rad[:2]
	gfun, hfun = [lambda x: ar_rad[2], lambda x: ar_rad[3]]
	_diff_spec = lambda y,x: diff_spec(*(args+(x,y)))*np.sin(x)
	return dblquad(_diff_spec, a, b, gfun, hfun, **kwargs)[0]
	
def int_spec_mat(diff_spec, omega = solid_angle(), *args, **kwargs):
	theta_inc = omega.theta_inc
	beta_inc = omega.beta_inc
	intmat = np.array([[int_spec(diff_spec,
			[theta-theta_inc/2.,theta+theta_inc/2.,beta-beta_inc/2.,beta+beta_inc/2.],
			*args, **kwargs) for theta in omega.theta_arr] for beta in omega.beta_arr])
	return intmat

	
def dcs_mat(dcs, omega = solid_angle(),*args):
	return dcs(*(args+(omega.theta_mesh, omega.beta_mesh)))

def dcs_mat_integral(dcs, omega = solid_angle(),*args):
	return dcs(*(args+(omega.theta_mesh, omega.beta_mesh)))*np.sin(omega.theta_mesh)*omega.delta_theta*omega.delta_beta

# dcs_mat_integral = int_spec_mat

def _spec_xrf_arr(energy_incident, Z, omega=None): # mac only
	"""
	returns the xrf spectrum from an ideal sample (MAC only) at solid angle omega
	"""
	if not omega:
		subtend = 4*np.pi
	else:
		subtend = omega.subtend
	# get xrf lines
	line = np.array(LineList2)
	energy = line_energy(Z,line)*1000
	# remove unavailable lines
	none_zero = np.where(energy>0)
	energy = energy[none_zero]
	line = line[none_zero]
	# get xrf line mac
	spec = mac_pe_line_kissel_cascade(energy_incident,Z,line)*subtend/(4*np.pi)
	# remove unavailable lines
	none_zero = np.where(spec>0)
	energy = energy[none_zero]
	spec = spec[none_zero]
	return energy, spec	

# # # compound, ideal sample, single layer
def spec_compton(compound, ev0, omega=solid_angle()):
	mac = dcs_mat_integral(compound.dmac_compton_pol, omega, ev0)
	energy_scattered = ev_scattered(ev0, omega.theta_rad_arr).reshape(1,-1)
	return energy_scattered, mac

def spec_rayleigh(compound, ev0, omega=solid_angle()):
	"""
	returns the Rayleigh spectrum from an ideal sample (dcs only) at solid angle omega
	"""
	mac = dcs_mat_integral(compound.dmac_rayleigh_pol, omega, ev0)
	ev0 = np.array(ev0,ndmin=2)
	return ev0, mac
	
def spec_xrf(compound, ev0, omega=solid_angle()):
	Z_arr = compound.Z_arr
	p_arr = compound.p_arr
	Z_out = []
	energy = []
	spec = []
	domega = np.sin(omega.theta_mesh)*omega.delta_theta*omega.delta_beta
	for i in range(len(Z_arr)):
		Z = Z_arr[i]
		p = p_arr[i]
		e, s = _spec_xrf_arr(ev0,Z,omega)
		if s.any():
			s *= p*atomic_weight(Z)/compound.molecular_weight
			e = e.reshape(-1,1,1)
			s = s.reshape(-1,1,1)*domega/(4*np.pi)
			energy.append(e)
			spec.append(s)
			Z_out.append(Z)
	return energy, spec, Z_out
		
# # # output in 'spectrum' format
# compound, ideal sample, single layer
def rawspec_xrf(compound, ev0, omega=solid_angle()):
	Z_arr = compound.Z_arr
	p_arr = compound.p_arr
	domega = np.sin(omega.theta_mesh)*omega.delta_theta*omega.delta_beta
	_out = []
	for i in range(len(Z_arr)):
		Z = Z_arr[i]
		p = p_arr[i]
		ev_arr, spec_arr = _spec_xrf_arr(ev0,Z,omega=None)
		if spec_arr.any():
			spec_arr *= p*atomic_weight(Z)/compound.molecular_weight
			ev_mat = ev_arr.reshape(-1,1,1)
			intensity_mat = spec_arr.reshape(-1,1,1)*domega/(4*np.pi)
			_out.append(spectrum(ev_mat = ev_mat, intensity_mat = intensity_mat,omega=deepcopy(omega),name=number2symbol(Z), stype = 'xrf'))
	return _out

	
def rawspec_compton(compound, ev0, omega=solid_angle()):
	mac_mat = dcs_mat_integral(compound.dmac_compton_pol, omega, ev0)
	energy_scattered_arr = ev_scattered(ev0, omega.theta_rad_arr)
	energy_scattered_mat = energy_scattered_arr.reshape(1,-1)
	return spectrum(ev_mat=energy_scattered_mat,intensity_mat=mac_mat,omega=deepcopy(omega),name='Compton', stype = 'comp')


def rawspec_rayleigh(compound, ev0, omega=solid_angle()):
	"""
	returns the Rayleigh spectrum from an ideal sample (dcs only) at solid angle omega
	"""
	mac_mat = dcs_mat_integral(compound.dmac_rayleigh_pol, omega, ev0)
	ev_mat = np.array(ev0,ndmin=2)
	return spectrum(ev_mat = ev_mat, intensity_mat = mac_mat, omega=deepcopy(omega), name='Rayleigh', stype = 'ray')

	
# # # sample geometry
def spec_monolayer(ev0, ev_mat, intensity_mat, psi_rad, psiprime_rad, mac_tot, density, thickness):
	temp = (mac_tot(ev0)/np.sin(psi_rad)+mac_tot(ev_mat)/np.sin(psiprime_rad))
	temp2 = mac_tot(ev_mat)/np.sin(psiprime_rad)
	refl = psiprime_rad > 0
	trans = psiprime_rad < 0
	refl_mat = 0
	trans_mat = 0
	if refl.any():
		refl_mat = refl*intensity_mat/(np.sin(psi_rad)*temp)*(1-np.exp(-temp*thickness*density))
	if trans.any():
		trans_mat = trans*intensity_mat/(np.sin(psi_rad)*temp)*(np.exp(temp2*thickness*density)-np.exp(-mac_tot(ev0)/np.sin(psi_rad)*density*thickness))
	intensity_mat = refl_mat + trans_mat
	intensity_mat[np.isnan(intensity_mat)] = 0
	return intensity_mat

def rawspec_multilayer(multilayer,illumination,omega=solid_angle()):
	ev0 = illumination.ev0
	psi_rad = illumination.psi_rad
	alpha_rad = illumination.alpha_rad
	psiprime_rad = psi_prime(omega.theta_mesh,omega.beta_mesh,psi_rad,alpha_rad)
	
	refl = psiprime_rad > 0
	trans = psiprime_rad < 0
	
	out = []
	xrf = []
	comp = []
	ray = []
	
	for sample in multilayer.sample_arr:
		# XRF spectrum
		xrf1 = rawspec_xrf(sample.compound, ev0, omega)
		# Rayleigh spectrum
		ray1 = rawspec_rayleigh(sample.compound, ev0, omega)
		# Compton spectrum
		comp1 = rawspec_compton(sample.compound, ev0, omega)
		# self-absorption within one layer
		for x in xrf1+[ray1,comp1]:
			intensity_mat = x.intensity_mat
			x.intensity_mat = spec_monolayer(ev0, x.ev_mat, x.intensity_mat, psi_rad, psiprime_rad, sample.compound.mac_total, sample.density, sample.thickness)
		# self-absorption from other layers
		for s1 in multilayer.sample_arr[:sample.layer]:
			for x in xrf1+[ray1,comp1]:
				x_refl = 0
				x_trans = 0
				if refl.any():
					x_refl = x.intensity_mat*refl*att_refl(ev0,x.ev_mat,psi_rad,psiprime_rad,s1.compound.mac_total,s1.density,s1.thickness)
				if trans.any():
					x_trans = x.intensity_mat*trans*att_trans_in(ev0,psi_rad,s1.compound.mac_total,s1.density,s1.thickness)
				x.intensity_mat = x_refl+x_trans
				x.intensity_mat[np.isnan(x.intensity_mat)] = 0
		if trans.any():
			for s2 in multilayer.sample_arr[sample.layer+1:]:
				for x in xrf1+[ray1,comp1]:
					x_trans = x.intensity_mat*trans*att_trans_out(x.ev_mat,psiprime_rad,s2.compound.mac_total,s2.density,s2.thickness)
					x.intensity_mat = x.intensity_mat*refl+x_trans
					x.intensity_mat[np.isnan(x.intensity_mat)] = 0
		xrf.append(xrf1)
		# Rayleigh spectrum
		ray.append(ray1)
		# Compton spectrum
		comp.append(comp1)
	_xrf = []
	_element_arr = []
	for x in flatten(xrf):
		try:
			 n = _element_arr.index(x.name)
			 _xrf[n].intensity_mat += x.intensity_mat
		except ValueError:
			_xrf.append(x)
			_element_arr.append(x.name)
	xrf = _xrf
	return xrf, ray, comp

def att_refl(ev0,ev,psi_rad,psiprime_rad,mac_tot,density,thickness):
	return np.exp(-(mac_tot(ev0)/np.sin(psi_rad)+mac_tot(ev)/np.sin(psiprime_rad))*density*thickness)

def att_trans_in(ev0,psi_rad,mac_tot,density,thickness):
	return np.exp(-mac_tot(ev0)/np.sin(psi_rad)*density*thickness)
	
def att_trans_out(ev,psiprime_rad,mac_tot,density,thickness):
	return np.exp(mac_tot(ev)/np.sin(psiprime_rad)*density*thickness)
	