import numpy as np
from scipy.special import erf
from copy import deepcopy

###constants
r_e = 2.8179403267e-13 #classical radius of electron in cm
m_e = 0.510998e6 #electron rest mass in eV/c^2
m_p = 938.272046e6 #proton rest mass in eV/c^2
a_H = 5.2917721092e-9 #Bohr radius in cm
hbar = 6.58211928e-16 #hbar in eV*s
c = 2.99792458e10	#speed of light in cm/s
alpha_f = hbar*c/(m_e*a_H) #fine structure constant
a_TF = (9*np.pi**2/2**7)**(1/3.)*a_H
N_A = 6.0221413e+23 #Avogadro constant

def erfc(x):
	return 1 - erf(x)

def unit(*args):
	return 1

def deg2rad(deg):
	"""
	Same as np.radians()
	"""
	deg = np.asarray(deg)
	return deg*np.pi/180.

def rad2deg(rad):
	"""
	Same as np.degrees()
	"""
	rad = np.asarray(rad)
	return rad/np.pi*180.

def ev2nm(ev):
	"""
	Convert photon energy (in ev) to wavelength (in nm).
	"""
	return 1239.84187/(1.*ev)
	
def theta2x(theta_rad,ev):
	"""
	theta in rad -> x in cm^-1
	"""
	return np.sin(theta_rad/2.)/(1.e-7*ev2nm(ev))


def arctan2(x):
	# like np.arctan, but the returned value is in [0,np.pi]
	rad = np.arctan(x)
	if rad<0:
		rad += np.pi
	return rad
	
def spherical2cartesian(r,theta_rad,beta_rad):
	"""
	Convert spherical coordinates (r,theta_rad,beta_rad) to Cartesian coordinates (x,y,z).
	"""
	x = r*np.sin(theta_rad)*np.cos(beta_rad)
	y = r*np.sin(theta_rad)*np.sin(beta_rad)
	z = r*np.cos(theta_rad)
	return x,y,z

def cartesian2spherical(x,y,z):
	"""
	Convert Cartesian coordinates (x,y,z) to spherical coordinates (r,theta_rad,beta_rad).
	"""
	r = np.sqrt(x**2+y**2+z**2)
	theta_rad = arctan2(np.sqrt(x**2+y**2)/(1.*z))
	beta_rad = arctan2(y/(1.*x))
	return r, theta_rad, beta_rad

def rotationy(rad):
	#return the rotation matrix for rotation about y axis
	return np.array([[np.cos(rad),0,np.sin(rad)],[0,1,0],[-np.sin(rad),0,np.cos(rad)]])

def rotationx(rad):
	#return the rotation matrix for rotation about x axis
	return np.array([[1,0,0],[0,np.cos(rad),-np.sin(rad)],[0,np.sin(rad),np.cos(rad)]])

	
def psi_prime_1(theta,beta,psi,alpha=0):
	# all angles in radians
	theta0 = np.pi/2.+psi
	theta_prime = cartesian2spherical(*np.dot(np.dot(np.transpose(rotationx(alpha)),np.transpose(rotationy(theta0))),np.array(spherical2cartesian(1,theta,beta))))[1]
	return np.pi/2.-theta_prime
	
def psi_prime(theta,beta,psi,alpha=0):
	return np.array(np.frompyfunc(psi_prime_1,4,1)(theta,beta,psi,alpha),dtype='float')
	
def flatten(lst):
	"""
	Flatten a list.
	"""
	return sum( ([x] if not isinstance(x, list) else flatten(x)
		     for x in lst), [] )

			 
def compound_parser(CP):
	Z_arr = []
	p_arr = []
	p = ''
	Z = CP[0]
	CP += 'A'
	for x in CP[1:]:
		if x.isdigit():
			p += x
		elif x.islower():
			Z += x
		elif x.isupper():
			Z_arr.append(Z)
			try:
				p_arr.append(int(p))
			except ValueError:
				p_arr.append(1)				
			p = ''
			Z = x
	return Z_arr, p_arr