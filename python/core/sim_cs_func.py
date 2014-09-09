# cs_func
from sim_common import *
from sim_xraylib import *

# # # Rayleigh scattering cross sections # # #

###unpolarized cross section

###polarized cross section

##differential cross section for single electron

def thomson_pol(theta_rad, beta_rad):
	return r_e**2*(np.cos(beta_rad)**2*np.cos(theta_rad)**2+np.sin(beta_rad)**2)

def thomson_unpol(theta_rad):
	return r_e**2*(1+np.cos(theta_rad)**2)/2.

def dcs_rayleigh_pol(ev, Z, theta_rad, beta_rad):
	dcs = aff_xraylib(1e-8*theta2x(theta_rad, ev), Z)**2*thomson_pol(theta_rad, beta_rad)
	return dcs
	
def dcs_rayleigh_unpol(ev, Z, theta_rad):
	dcs = aff_xraylib(1e-8*theta2x(theta_rad, ev), Z)**2*thomson_unpol(theta_rad)
	return dcs
	

	
# # # Compton scattering cross sections # # #
	
###compton dcs of single electron
##dcs
def ev_scattered(ev, theta_rad):
	"""
	return the energy in ev of the Compton scattered photons.
	"""
	return ev/(1+ev/m_e*(1-np.cos(theta_rad)))

def klein_nishina_unpol(ev,theta_rad):
	"""
	Differential Klein-Nishina cross section in cm^2 for unpolarized radiation.
	"""
	return r_e**2/2.*(ev_scattered(ev, theta_rad)/ev)**2*(ev_scattered(ev, theta_rad)/ev+ev/ev_scattered(ev, theta_rad)-np.sin(theta_rad)**2)

def klein_nishina_pol(ev, theta_rad, beta_rad):
	"""
	Differential Klein-Nishina cross section in cm^2 for polarized radiation.
	"""
	return r_e**2/2.*(ev_scattered(ev, theta_rad)/ev)**2*(ev_scattered(ev, theta_rad)/ev+ev/ev_scattered(ev, theta_rad)-2*np.sin(theta_rad)**2*np.cos(beta_rad)**2)

##total cs
def klein_nishina_total_col(ev):
	"""
	Total Klein-Nishina collision cross section in cm^2.
	"""
	a = ev/m_e
	return 2*np.pi*r_e**2*((1+a)/(a**3)*(2*a*(1+a)/(1+2*a)-np.log(1+2*a))+np.log(1+2*a)/(2*a)-(1+3*a)/(1+2*a)**2)

def klein_nishina_total_sca(ev):
	"""
	Total Klein-Nishina scattering cross section in cm^2.
	"""
	a = ev/m_e
	return np.pi*r_e**2*(np.log(1+2*a)/(a**3)+2*(1+a)*(2*a**2-2*a-1)/(a**2*(1+2*a)**2)+8*a**2/(3*(1+2*a)**3))

###compton dcs of atom
def dcs_compton_pol(ev, Z, theta_rad, beta_rad):
	"""
	Differential Compton cross section in cm^2 for polarized radiation.
	"""
	return klein_nishina_pol(ev, theta_rad, beta_rad)*isf_xraylib(theta2x(theta_rad,ev)*1e-8,Z)
	
def dcs_compton_unpol(ev, Z, theta_rad):
	"""
	Differential Compton cross section in cm^2 for unpolarized radiation.
	"""
	return klein_nishina_unpol(ev, theta_rad)*isf_xraylib(theta2x(theta_rad,ev)*1e-8,Z)

