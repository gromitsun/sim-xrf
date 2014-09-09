import numpy as np
from xraylib import *


ray_dcs_pol=np.frompyfunc(DCSPb_Rayl,4,1)


def aff_xraylib(x,Z):
	"""
	Atomic form factor.
	x is in angstrom^-1
	"""
	return np.array(np.frompyfunc(FF_Rayl,2,1)(Z,x),dtype='float')

def isf_xraylib(x,Z):
	"""
	Incoherent scattering function.
	x is in angstrom^-1
	"""
	return np.array(np.frompyfunc(SF_Compt,2,1)(Z,x),dtype='float')


# cross sections	

def mac_pe_line(ev0, Z, line):
	return np.array(np.frompyfunc(CS_FluorLine,3,1)(Z, line, ev0/1000.),dtype='float')
	
def mac_pe_line_kissel_cascade(ev0, Z, line):
	return np.array(np.frompyfunc(CS_FluorLine_Kissel_Cascade,3,1)(Z, line, ev0/1000.),dtype='float')

def mac_pe_line_kissel(ev0, Z, line):
	return np.array(np.frompyfunc(CS_FluorLine_Kissel,3,1)(Z, line, ev0/1000.),dtype='float')	

	
def mac_pe(ev0, Z):
	return np.array(np.frompyfunc(CS_Photo,2,1)(Z, ev0/1000.),dtype='float')

def mac_total(ev0, Z):
	return np.array(np.frompyfunc(CS_Total,2,1)(Z, ev0/1000.),dtype='float')

def mac_compton(ev0, Z):
	return np.array(np.frompyfunc(CS_Compt,2,1)(Z, ev0/1000.),dtype='float')
	
def mac_rayleigh(ev0, Z):
	return np.array(np.frompyfunc(CS_Rayl,2,1)(Z, ev0/1000.),dtype='float')


def atomic_weight(Z):
	return np.array(np.frompyfunc(AtomicWeight,1,1)(Z),dtype='float')


def molecular_weight(Z_arr,p_arr):
	M = 0
	for i in range(len(Z_arr)):
		Z = symbol2number(Z_arr[i])
		p = p_arr[i]
		M += atomic_weight(Z)*p
	return M
	
def element_density(Z):
	return np.array(np.frompyfunc(ElementDensity,1,1)(Z),dtype='float')
	
def xrf_yield(Z, shell):
	return np.array(np.frompyfunc(FluorYield,2,1)(Z,shell), dtype = 'float')

def line_energy(Z, line):
	return np.array(np.frompyfunc(LineEnergy,2,1)(Z,line), dtype = 'float')


	
# compound processing
def number2symbol(Z):
	try:
		return AtomicNumberToSymbol(int(float(Z)))
	except ValueError:
		return Z
	except TypeError:
		out = []
		for x in Z:
			out.append(number2symbol(x))
		return out

def symbol2number(Z):
	try:
		return int(Z)
	except ValueError:
		return SymbolToAtomicNumber(Z)
	except TypeError:
		out = []
		for x in Z:
			out.append(symbol2number(x))
		return out



KA_LINE = 0 
KB_LINE = 1 
LA_LINE = 2 
LB_LINE = 3 
LG_LINE = 4 
LB2_LINE = 5
MA_LINE = 6 


LineList0=[
KA_LINE, # 0 
KB_LINE, # 1 
LA_LINE, # 2 
LB_LINE, # 3 
LG_LINE, # 4 
# LB2_LINE, # 5
MA_LINE # 6 
]

LineList1=[
KA1_LINE, #_xraylib.KA1_LINE,
KA2_LINE, #_xraylib.KA2_LINE,
KB1_LINE, #_xraylib.KB1_LINE,
KB2_LINE, #_xraylib.KB2_LINE,
KB3_LINE, #_xraylib.KB3_LINE,
KB4_LINE, #_xraylib.KB4_LINE,
KB5_LINE, #_xraylib.KB5_LINE,
LA1_LINE, #_xraylib.LA1_LINE,
LA2_LINE, #_xraylib.LA2_LINE,
LB1_LINE, #_xraylib.LB1_LINE,
LB2_LINE, #_xraylib.LB2_LINE,
LB3_LINE, #_xraylib.LB3_LINE,
LB4_LINE, #_xraylib.LB4_LINE,
LB5_LINE, #_xraylib.LB5_LINE,
LB6_LINE, #_xraylib.LB6_LINE,
LB7_LINE, #_xraylib.LB7_LINE,
LB9_LINE, #_xraylib.LB9_LINE,
LB10_LINE, #_xraylib.LB10_LINE,
LB15_LINE, #_xraylib.LB15_LINE,
LB17_LINE, #_xraylib.LB17_LINE,
LG1_LINE, #_xraylib.LG1_LINE,
LG2_LINE, #_xraylib.LG2_LINE,
LG3_LINE, #_xraylib.LG3_LINE,
LG4_LINE, #_xraylib.LG4_LINE,
LG5_LINE, #_xraylib.LG5_LINE,
LG6_LINE, #_xraylib.LG6_LINE,
LG8_LINE, #_xraylib.LG8_LINE,
LE_LINE, #_xraylib.LE_LINE,
LL_LINE, #_xraylib.LL_LINE,
LS_LINE, #_xraylib.LS_LINE,
LT_LINE, #_xraylib.LT_LINE,
LU_LINE, #_xraylib.LU_LINE,
LV_LINE, #_xraylib.LV_LINE,
MA1_LINE, #_xraylib.MA1_LINE,
MA2_LINE, #_xraylib.MA2_LINE,
MB_LINE, #_xraylib.MB_LINE,
MG_LINE] #_xraylib.MG_LINE,

LineList2=[
KL1_LINE, # _xraylib.KL1_LINE
KL2_LINE, # _xraylib.KL2_LINE
KL3_LINE, # _xraylib.KL3_LINE
KM1_LINE, # _xraylib.KM1_LINE
KM2_LINE, # _xraylib.KM2_LINE
KM3_LINE, # _xraylib.KM3_LINE
KM4_LINE, # _xraylib.KM4_LINE
KM5_LINE, # _xraylib.KM5_LINE
KN1_LINE, # _xraylib.KN1_LINE
KN2_LINE, # _xraylib.KN2_LINE
KN3_LINE, # _xraylib.KN3_LINE
KN4_LINE, # _xraylib.KN4_LINE
KN5_LINE, # _xraylib.KN5_LINE
KN6_LINE, # _xraylib.KN6_LINE
KN7_LINE, # _xraylib.KN7_LINE
KO_LINE, # _xraylib.KO_LINE
KO1_LINE, # _xraylib.KO1_LINE
KO2_LINE, # _xraylib.KO2_LINE
KO3_LINE, # _xraylib.KO3_LINE
KO4_LINE, # _xraylib.KO4_LINE
KO5_LINE, # _xraylib.KO5_LINE
KO6_LINE, # _xraylib.KO6_LINE
KO7_LINE, # _xraylib.KO7_LINE
KP_LINE, # _xraylib.KP_LINE
KP1_LINE, # _xraylib.KP1_LINE
KP2_LINE, # _xraylib.KP2_LINE
KP3_LINE, # _xraylib.KP3_LINE
KP4_LINE, # _xraylib.KP4_LINE
KP5_LINE, # _xraylib.KP5_LINE
L1L2_LINE, # _xraylib.L1L2_LINE
L1L3_LINE, # _xraylib.L1L3_LINE
L1M1_LINE, # _xraylib.L1M1_LINE
L1M2_LINE, # _xraylib.L1M2_LINE
L1M3_LINE, # _xraylib.L1M3_LINE
L1M4_LINE, # _xraylib.L1M4_LINE
L1M5_LINE, # _xraylib.L1M5_LINE
L1N1_LINE, # _xraylib.L1N1_LINE
L1N2_LINE, # _xraylib.L1N2_LINE
L1N3_LINE, # _xraylib.L1N3_LINE
L1N4_LINE, # _xraylib.L1N4_LINE
L1N5_LINE, # _xraylib.L1N5_LINE
L1N6_LINE, # _xraylib.L1N6_LINE
L1N67_LINE, # _xraylib.L1N67_LINE
L1N7_LINE, # _xraylib.L1N7_LINE
L1O1_LINE, # _xraylib.L1O1_LINE
L1O2_LINE, # _xraylib.L1O2_LINE
L1O3_LINE, # _xraylib.L1O3_LINE
L1O4_LINE, # _xraylib.L1O4_LINE
L1O45_LINE, # _xraylib.L1O45_LINE
L1O5_LINE, # _xraylib.L1O5_LINE
L1O6_LINE, # _xraylib.L1O6_LINE
L1O7_LINE, # _xraylib.L1O7_LINE
L1P1_LINE, # _xraylib.L1P1_LINE
L1P2_LINE, # _xraylib.L1P2_LINE
L1P23_LINE, # _xraylib.L1P23_LINE
L1P3_LINE, # _xraylib.L1P3_LINE
L1P4_LINE, # _xraylib.L1P4_LINE
L1P5_LINE, # _xraylib.L1P5_LINE
L2L3_LINE, # _xraylib.L2L3_LINE
L2M1_LINE, # _xraylib.L2M1_LINE
L2M2_LINE, # _xraylib.L2M2_LINE
L2M3_LINE, # _xraylib.L2M3_LINE
L2M4_LINE, # _xraylib.L2M4_LINE
L2M5_LINE, # _xraylib.L2M5_LINE
L2N1_LINE, # _xraylib.L2N1_LINE
L2N2_LINE, # _xraylib.L2N2_LINE
L2N3_LINE, # _xraylib.L2N3_LINE
L2N4_LINE, # _xraylib.L2N4_LINE
L2N5_LINE, # _xraylib.L2N5_LINE
L2N6_LINE, # _xraylib.L2N6_LINE
L2N7_LINE, # _xraylib.L2N7_LINE
L2O1_LINE, # _xraylib.L2O1_LINE
L2O2_LINE, # _xraylib.L2O2_LINE
L2O3_LINE, # _xraylib.L2O3_LINE
L2O4_LINE, # _xraylib.L2O4_LINE
L2O5_LINE, # _xraylib.L2O5_LINE
L2O6_LINE, # _xraylib.L2O6_LINE
L2O7_LINE, # _xraylib.L2O7_LINE
L2P1_LINE, # _xraylib.L2P1_LINE
L2P2_LINE, # _xraylib.L2P2_LINE
L2P23_LINE, # _xraylib.L2P23_LINE
L2P3_LINE, # _xraylib.L2P3_LINE
L2P4_LINE, # _xraylib.L2P4_LINE
L2P5_LINE, # _xraylib.L2P5_LINE
L2Q1_LINE, # _xraylib.L2Q1_LINE
L3M1_LINE, # _xraylib.L3M1_LINE
L3M2_LINE, # _xraylib.L3M2_LINE
L3M3_LINE, # _xraylib.L3M3_LINE
L3M4_LINE, # _xraylib.L3M4_LINE
L3M5_LINE, # _xraylib.L3M5_LINE
L3N1_LINE, # _xraylib.L3N1_LINE
L3N2_LINE, # _xraylib.L3N2_LINE
L3N3_LINE, # _xraylib.L3N3_LINE
L3N4_LINE, # _xraylib.L3N4_LINE
L3N5_LINE, # _xraylib.L3N5_LINE
L3N6_LINE, # _xraylib.L3N6_LINE
L3N7_LINE, # _xraylib.L3N7_LINE
L3O1_LINE, # _xraylib.L3O1_LINE
L3O2_LINE, # _xraylib.L3O2_LINE
L3O3_LINE, # _xraylib.L3O3_LINE
L3O4_LINE, # _xraylib.L3O4_LINE
L3O45_LINE, # _xraylib.L3O45_LINE
L3O5_LINE, # _xraylib.L3O5_LINE
L3O6_LINE, # _xraylib.L3O6_LINE
L3O7_LINE, # _xraylib.L3O7_LINE
L3P1_LINE, # _xraylib.L3P1_LINE
L3P2_LINE, # _xraylib.L3P2_LINE
L3P23_LINE, # _xraylib.L3P23_LINE
L3P3_LINE, # _xraylib.L3P3_LINE
L3P4_LINE, # _xraylib.L3P4_LINE
L3P45_LINE, # _xraylib.L3P45_LINE
L3P5_LINE, # _xraylib.L3P5_LINE
L3Q1_LINE, # _xraylib.L3Q1_LINE
M1M2_LINE, # _xraylib.M1M2_LINE
M1M3_LINE, # _xraylib.M1M3_LINE
M1M4_LINE, # _xraylib.M1M4_LINE
M1M5_LINE, # _xraylib.M1M5_LINE
M1N1_LINE, # _xraylib.M1N1_LINE
M1N2_LINE, # _xraylib.M1N2_LINE
M1N3_LINE, # _xraylib.M1N3_LINE
M1N4_LINE, # _xraylib.M1N4_LINE
M1N5_LINE, # _xraylib.M1N5_LINE
M1N6_LINE, # _xraylib.M1N6_LINE
M1N7_LINE, # _xraylib.M1N7_LINE
M1O1_LINE, # _xraylib.M1O1_LINE
M1O2_LINE, # _xraylib.M1O2_LINE
M1O3_LINE, # _xraylib.M1O3_LINE
M1O4_LINE, # _xraylib.M1O4_LINE
M1O5_LINE, # _xraylib.M1O5_LINE
M1O6_LINE, # _xraylib.M1O6_LINE
M1O7_LINE, # _xraylib.M1O7_LINE
M1P1_LINE, # _xraylib.M1P1_LINE
M1P2_LINE, # _xraylib.M1P2_LINE
M1P3_LINE, # _xraylib.M1P3_LINE
M1P4_LINE, # _xraylib.M1P4_LINE
M1P5_LINE, # _xraylib.M1P5_LINE
M2M3_LINE, # _xraylib.M2M3_LINE
M2M4_LINE, # _xraylib.M2M4_LINE
M2M5_LINE, # _xraylib.M2M5_LINE
M2N1_LINE, # _xraylib.M2N1_LINE
M2N2_LINE, # _xraylib.M2N2_LINE
M2N3_LINE, # _xraylib.M2N3_LINE
M2N4_LINE, # _xraylib.M2N4_LINE
M2N5_LINE, # _xraylib.M2N5_LINE
M2N6_LINE, # _xraylib.M2N6_LINE
M2N7_LINE, # _xraylib.M2N7_LINE
M2O1_LINE, # _xraylib.M2O1_LINE
M2O2_LINE, # _xraylib.M2O2_LINE
M2O3_LINE, # _xraylib.M2O3_LINE
M2O4_LINE, # _xraylib.M2O4_LINE
M2O5_LINE, # _xraylib.M2O5_LINE
M2O6_LINE, # _xraylib.M2O6_LINE
M2O7_LINE, # _xraylib.M2O7_LINE
M2P1_LINE, # _xraylib.M2P1_LINE
M2P2_LINE, # _xraylib.M2P2_LINE
M2P3_LINE, # _xraylib.M2P3_LINE
M2P4_LINE, # _xraylib.M2P4_LINE
M2P5_LINE, # _xraylib.M2P5_LINE
M3M4_LINE, # _xraylib.M3M4_LINE
M3M5_LINE, # _xraylib.M3M5_LINE
M3N1_LINE, # _xraylib.M3N1_LINE
M3N2_LINE, # _xraylib.M3N2_LINE
M3N3_LINE, # _xraylib.M3N3_LINE
M3N4_LINE, # _xraylib.M3N4_LINE
M3N5_LINE, # _xraylib.M3N5_LINE
M3N6_LINE, # _xraylib.M3N6_LINE
M3N7_LINE, # _xraylib.M3N7_LINE
M3O1_LINE, # _xraylib.M3O1_LINE
M3O2_LINE, # _xraylib.M3O2_LINE
M3O3_LINE, # _xraylib.M3O3_LINE
M3O4_LINE, # _xraylib.M3O4_LINE
M3O5_LINE, # _xraylib.M3O5_LINE
M3O6_LINE, # _xraylib.M3O6_LINE
M3O7_LINE, # _xraylib.M3O7_LINE
M3P1_LINE, # _xraylib.M3P1_LINE
M3P2_LINE, # _xraylib.M3P2_LINE
M3P3_LINE, # _xraylib.M3P3_LINE
M3P4_LINE, # _xraylib.M3P4_LINE
M3P5_LINE, # _xraylib.M3P5_LINE
M3Q1_LINE, # _xraylib.M3Q1_LINE
M4M5_LINE, # _xraylib.M4M5_LINE
M4N1_LINE, # _xraylib.M4N1_LINE
M4N2_LINE, # _xraylib.M4N2_LINE
M4N3_LINE, # _xraylib.M4N3_LINE
M4N4_LINE, # _xraylib.M4N4_LINE
M4N5_LINE, # _xraylib.M4N5_LINE
M4N6_LINE, # _xraylib.M4N6_LINE
M4N7_LINE, # _xraylib.M4N7_LINE
M4O1_LINE, # _xraylib.M4O1_LINE
M4O2_LINE, # _xraylib.M4O2_LINE
M4O3_LINE, # _xraylib.M4O3_LINE
M4O4_LINE, # _xraylib.M4O4_LINE
M4O5_LINE, # _xraylib.M4O5_LINE
M4O6_LINE, # _xraylib.M4O6_LINE
M4O7_LINE, # _xraylib.M4O7_LINE
M4P1_LINE, # _xraylib.M4P1_LINE
M4P2_LINE, # _xraylib.M4P2_LINE
M4P3_LINE, # _xraylib.M4P3_LINE
M4P4_LINE, # _xraylib.M4P4_LINE
M4P5_LINE, # _xraylib.M4P5_LINE
M5N1_LINE, # _xraylib.M5N1_LINE
M5N2_LINE, # _xraylib.M5N2_LINE
M5N3_LINE, # _xraylib.M5N3_LINE
M5N4_LINE, # _xraylib.M5N4_LINE
M5N5_LINE, # _xraylib.M5N5_LINE
M5N6_LINE, # _xraylib.M5N6_LINE
M5N7_LINE, # _xraylib.M5N7_LINE
M5O1_LINE, # _xraylib.M5O1_LINE
M5O2_LINE, # _xraylib.M5O2_LINE
M5O3_LINE, # _xraylib.M5O3_LINE
M5O4_LINE, # _xraylib.M5O4_LINE
M5O5_LINE, # _xraylib.M5O5_LINE
M5O6_LINE, # _xraylib.M5O6_LINE
M5O7_LINE, # _xraylib.M5O7_LINE
M5P1_LINE, # _xraylib.M5P1_LINE
M5P2_LINE, # _xraylib.M5P2_LINE
M5P3_LINE, # _xraylib.M5P3_LINE
M5P4_LINE, # _xraylib.M5P4_LINE
M5P5_LINE, # _xraylib.M5P5_LINE
N1N2_LINE, # _xraylib.N1N2_LINE
N1N3_LINE, # _xraylib.N1N3_LINE
N1N4_LINE, # _xraylib.N1N4_LINE
N1N5_LINE, # _xraylib.N1N5_LINE
N1N6_LINE, # _xraylib.N1N6_LINE
N1N7_LINE, # _xraylib.N1N7_LINE
N1O1_LINE, # _xraylib.N1O1_LINE
N1O2_LINE, # _xraylib.N1O2_LINE
N1O3_LINE, # _xraylib.N1O3_LINE
N1O4_LINE, # _xraylib.N1O4_LINE
N1O5_LINE, # _xraylib.N1O5_LINE
N1O6_LINE, # _xraylib.N1O6_LINE
N1O7_LINE, # _xraylib.N1O7_LINE
N1P1_LINE, # _xraylib.N1P1_LINE
N1P2_LINE, # _xraylib.N1P2_LINE
N1P3_LINE, # _xraylib.N1P3_LINE
N1P4_LINE, # _xraylib.N1P4_LINE
N1P5_LINE, # _xraylib.N1P5_LINE
N2N3_LINE, # _xraylib.N2N3_LINE
N2N4_LINE, # _xraylib.N2N4_LINE
N2N5_LINE, # _xraylib.N2N5_LINE
N2N6_LINE, # _xraylib.N2N6_LINE
N2N7_LINE, # _xraylib.N2N7_LINE
N2O1_LINE, # _xraylib.N2O1_LINE
N2O2_LINE, # _xraylib.N2O2_LINE
N2O3_LINE, # _xraylib.N2O3_LINE
N2O4_LINE, # _xraylib.N2O4_LINE
N2O5_LINE, # _xraylib.N2O5_LINE
N2O6_LINE, # _xraylib.N2O6_LINE
N2O7_LINE, # _xraylib.N2O7_LINE
N2P1_LINE, # _xraylib.N2P1_LINE
N2P2_LINE, # _xraylib.N2P2_LINE
N2P3_LINE, # _xraylib.N2P3_LINE
N2P4_LINE, # _xraylib.N2P4_LINE
N2P5_LINE, # _xraylib.N2P5_LINE
N3N4_LINE, # _xraylib.N3N4_LINE
N3N5_LINE, # _xraylib.N3N5_LINE
N3N6_LINE, # _xraylib.N3N6_LINE
N3N7_LINE, # _xraylib.N3N7_LINE
N3O1_LINE, # _xraylib.N3O1_LINE
N3O2_LINE, # _xraylib.N3O2_LINE
N3O3_LINE, # _xraylib.N3O3_LINE
N3O4_LINE, # _xraylib.N3O4_LINE
N3O5_LINE, # _xraylib.N3O5_LINE
N3O6_LINE, # _xraylib.N3O6_LINE
N3O7_LINE, # _xraylib.N3O7_LINE
N3P1_LINE, # _xraylib.N3P1_LINE
N3P2_LINE, # _xraylib.N3P2_LINE
N3P3_LINE, # _xraylib.N3P3_LINE
N3P4_LINE, # _xraylib.N3P4_LINE
N3P5_LINE, # _xraylib.N3P5_LINE
N4N5_LINE, # _xraylib.N4N5_LINE
N4N6_LINE, # _xraylib.N4N6_LINE
N4N7_LINE, # _xraylib.N4N7_LINE
N4O1_LINE, # _xraylib.N4O1_LINE
N4O2_LINE, # _xraylib.N4O2_LINE
N4O3_LINE, # _xraylib.N4O3_LINE
N4O4_LINE, # _xraylib.N4O4_LINE
N4O5_LINE, # _xraylib.N4O5_LINE
N4O6_LINE, # _xraylib.N4O6_LINE
N4O7_LINE, # _xraylib.N4O7_LINE
N4P1_LINE, # _xraylib.N4P1_LINE
N4P2_LINE, # _xraylib.N4P2_LINE
N4P3_LINE, # _xraylib.N4P3_LINE
N4P4_LINE, # _xraylib.N4P4_LINE
N4P5_LINE, # _xraylib.N4P5_LINE
N5N6_LINE, # _xraylib.N5N6_LINE
N5N7_LINE, # _xraylib.N5N7_LINE
N5O1_LINE, # _xraylib.N5O1_LINE
N5O2_LINE, # _xraylib.N5O2_LINE
N5O3_LINE, # _xraylib.N5O3_LINE
N5O4_LINE, # _xraylib.N5O4_LINE
N5O5_LINE, # _xraylib.N5O5_LINE
N5O6_LINE, # _xraylib.N5O6_LINE
N5O7_LINE, # _xraylib.N5O7_LINE
N5P1_LINE, # _xraylib.N5P1_LINE
N5P2_LINE, # _xraylib.N5P2_LINE
N5P3_LINE, # _xraylib.N5P3_LINE
N5P4_LINE, # _xraylib.N5P4_LINE
N5P5_LINE, # _xraylib.N5P5_LINE
N6N7_LINE, # _xraylib.N6N7_LINE
N6O1_LINE, # _xraylib.N6O1_LINE
N6O2_LINE, # _xraylib.N6O2_LINE
N6O3_LINE, # _xraylib.N6O3_LINE
N6O4_LINE, # _xraylib.N6O4_LINE
N6O5_LINE, # _xraylib.N6O5_LINE
N6O6_LINE, # _xraylib.N6O6_LINE
N6O7_LINE, # _xraylib.N6O7_LINE
N6P1_LINE, # _xraylib.N6P1_LINE
N6P2_LINE, # _xraylib.N6P2_LINE
N6P3_LINE, # _xraylib.N6P3_LINE
N6P4_LINE, # _xraylib.N6P4_LINE
N6P5_LINE, # _xraylib.N6P5_LINE
N7O1_LINE, # _xraylib.N7O1_LINE
N7O2_LINE, # _xraylib.N7O2_LINE
N7O3_LINE, # _xraylib.N7O3_LINE
N7O4_LINE, # _xraylib.N7O4_LINE
N7O5_LINE, # _xraylib.N7O5_LINE
N7O6_LINE, # _xraylib.N7O6_LINE
N7O7_LINE, # _xraylib.N7O7_LINE
N7P1_LINE, # _xraylib.N7P1_LINE
N7P2_LINE, # _xraylib.N7P2_LINE
N7P3_LINE, # _xraylib.N7P3_LINE
N7P4_LINE, # _xraylib.N7P4_LINE
N7P5_LINE, # _xraylib.N7P5_LINE
O1O2_LINE, # _xraylib.O1O2_LINE
O1O3_LINE, # _xraylib.O1O3_LINE
O1O4_LINE, # _xraylib.O1O4_LINE
O1O5_LINE, # _xraylib.O1O5_LINE
O1O6_LINE, # _xraylib.O1O6_LINE
O1O7_LINE, # _xraylib.O1O7_LINE
O1P1_LINE, # _xraylib.O1P1_LINE
O1P2_LINE, # _xraylib.O1P2_LINE
O1P3_LINE, # _xraylib.O1P3_LINE
O1P4_LINE, # _xraylib.O1P4_LINE
O1P5_LINE, # _xraylib.O1P5_LINE
O2O3_LINE, # _xraylib.O2O3_LINE
O2O4_LINE, # _xraylib.O2O4_LINE
O2O5_LINE, # _xraylib.O2O5_LINE
O2O6_LINE, # _xraylib.O2O6_LINE
O2O7_LINE, # _xraylib.O2O7_LINE
O2P1_LINE, # _xraylib.O2P1_LINE
O2P2_LINE, # _xraylib.O2P2_LINE
O2P3_LINE, # _xraylib.O2P3_LINE
O2P4_LINE, # _xraylib.O2P4_LINE
O2P5_LINE, # _xraylib.O2P5_LINE
O3O4_LINE, # _xraylib.O3O4_LINE
O3O5_LINE, # _xraylib.O3O5_LINE
O3O6_LINE, # _xraylib.O3O6_LINE
O3O7_LINE, # _xraylib.O3O7_LINE
O3P1_LINE, # _xraylib.O3P1_LINE
O3P2_LINE, # _xraylib.O3P2_LINE
O3P3_LINE, # _xraylib.O3P3_LINE
O3P4_LINE, # _xraylib.O3P4_LINE
O3P5_LINE, # _xraylib.O3P5_LINE
O4O5_LINE, # _xraylib.O4O5_LINE
O4O6_LINE, # _xraylib.O4O6_LINE
O4O7_LINE, # _xraylib.O4O7_LINE
O4P1_LINE, # _xraylib.O4P1_LINE
O4P2_LINE, # _xraylib.O4P2_LINE
O4P3_LINE, # _xraylib.O4P3_LINE
O4P4_LINE, # _xraylib.O4P4_LINE
O4P5_LINE, # _xraylib.O4P5_LINE
O5O6_LINE, # _xraylib.O5O6_LINE
O5O7_LINE, # _xraylib.O5O7_LINE
O5P1_LINE, # _xraylib.O5P1_LINE
O5P2_LINE, # _xraylib.O5P2_LINE
O5P3_LINE, # _xraylib.O5P3_LINE
O5P4_LINE, # _xraylib.O5P4_LINE
O5P5_LINE, # _xraylib.O5P5_LINE
O6O7_LINE, # _xraylib.O6O7_LINE
O6P4_LINE, # _xraylib.O6P4_LINE
O6P5_LINE, # _xraylib.O6P5_LINE
O7P4_LINE, # _xraylib.O7P4_LINE
O7P5_LINE, # _xraylib.O7P5_LINE
P1P2_LINE, # _xraylib.P1P2_LINE
P1P3_LINE, # _xraylib.P1P3_LINE
P1P4_LINE, # _xraylib.P1P4_LINE
P1P5_LINE, # _xraylib.P1P5_LINE
P2P3_LINE, # _xraylib.P2P3_LINE
P2P4_LINE, # _xraylib.P2P4_LINE
P2P5_LINE, # _xraylib.P2P5_LINE
P3P4_LINE, # _xraylib.P3P4_LINE
P3P5_LINE, # _xraylib.P3P5_LINE
]

ShellList = [
K_SHELL, #  0
L1_SHELL, # 1
L2_SHELL, # 2
L3_SHELL, # 3
M1_SHELL, # 4
M2_SHELL, # 5
M3_SHELL, # 6
M4_SHELL, # 7
M5_SHELL # 8 
]



K_LINES = [
KL1_LINE, # _xraylib.KL1_LINE
KL2_LINE, # _xraylib.KL2_LINE
KL3_LINE, # _xraylib.KL3_LINE
KM1_LINE, # _xraylib.KM1_LINE
KM2_LINE, # _xraylib.KM2_LINE
KM3_LINE, # _xraylib.KM3_LINE
KM4_LINE, # _xraylib.KM4_LINE
KM5_LINE, # _xraylib.KM5_LINE
KN1_LINE, # _xraylib.KN1_LINE
KN2_LINE, # _xraylib.KN2_LINE
KN3_LINE, # _xraylib.KN3_LINE
KN4_LINE, # _xraylib.KN4_LINE
KN5_LINE, # _xraylib.KN5_LINE
KN6_LINE, # _xraylib.KN6_LINE
KN7_LINE, # _xraylib.KN7_LINE
KO_LINE, # _xraylib.KO_LINE
KO1_LINE, # _xraylib.KO1_LINE
KO2_LINE, # _xraylib.KO2_LINE
KO3_LINE, # _xraylib.KO3_LINE
KO4_LINE, # _xraylib.KO4_LINE
KO5_LINE, # _xraylib.KO5_LINE
KO6_LINE, # _xraylib.KO6_LINE
KO7_LINE, # _xraylib.KO7_LINE
KP_LINE, # _xraylib.KP_LINE
KP1_LINE, # _xraylib.KP1_LINE
KP2_LINE, # _xraylib.KP2_LINE
KP3_LINE, # _xraylib.KP3_LINE
KP4_LINE, # _xraylib.KP4_LINE
KP5_LINE]

L_LINES = [
L1L2_LINE, # _xraylib.L1L2_LINE
L1L3_LINE, # _xraylib.L1L3_LINE
L1M1_LINE, # _xraylib.L1M1_LINE
L1M2_LINE, # _xraylib.L1M2_LINE
L1M3_LINE, # _xraylib.L1M3_LINE
L1M4_LINE, # _xraylib.L1M4_LINE
L1M5_LINE, # _xraylib.L1M5_LINE
L1N1_LINE, # _xraylib.L1N1_LINE
L1N2_LINE, # _xraylib.L1N2_LINE
L1N3_LINE, # _xraylib.L1N3_LINE
L1N4_LINE, # _xraylib.L1N4_LINE
L1N5_LINE, # _xraylib.L1N5_LINE
L1N6_LINE, # _xraylib.L1N6_LINE
L1N67_LINE, # _xraylib.L1N67_LINE
L1N7_LINE, # _xraylib.L1N7_LINE
L1O1_LINE, # _xraylib.L1O1_LINE
L1O2_LINE, # _xraylib.L1O2_LINE
L1O3_LINE, # _xraylib.L1O3_LINE
L1O4_LINE, # _xraylib.L1O4_LINE
L1O45_LINE, # _xraylib.L1O45_LINE
L1O5_LINE, # _xraylib.L1O5_LINE
L1O6_LINE, # _xraylib.L1O6_LINE
L1O7_LINE, # _xraylib.L1O7_LINE
L1P1_LINE, # _xraylib.L1P1_LINE
L1P2_LINE, # _xraylib.L1P2_LINE
L1P23_LINE, # _xraylib.L1P23_LINE
L1P3_LINE, # _xraylib.L1P3_LINE
L1P4_LINE, # _xraylib.L1P4_LINE
L1P5_LINE, # _xraylib.L1P5_LINE
L2L3_LINE, # _xraylib.L2L3_LINE
L2M1_LINE, # _xraylib.L2M1_LINE
L2M2_LINE, # _xraylib.L2M2_LINE
L2M3_LINE, # _xraylib.L2M3_LINE
L2M4_LINE, # _xraylib.L2M4_LINE
L2M5_LINE, # _xraylib.L2M5_LINE
L2N1_LINE, # _xraylib.L2N1_LINE
L2N2_LINE, # _xraylib.L2N2_LINE
L2N3_LINE, # _xraylib.L2N3_LINE
L2N4_LINE, # _xraylib.L2N4_LINE
L2N5_LINE, # _xraylib.L2N5_LINE
L2N6_LINE, # _xraylib.L2N6_LINE
L2N7_LINE, # _xraylib.L2N7_LINE
L2O1_LINE, # _xraylib.L2O1_LINE
L2O2_LINE, # _xraylib.L2O2_LINE
L2O3_LINE, # _xraylib.L2O3_LINE
L2O4_LINE, # _xraylib.L2O4_LINE
L2O5_LINE, # _xraylib.L2O5_LINE
L2O6_LINE, # _xraylib.L2O6_LINE
L2O7_LINE, # _xraylib.L2O7_LINE
L2P1_LINE, # _xraylib.L2P1_LINE
L2P2_LINE, # _xraylib.L2P2_LINE
L2P23_LINE, # _xraylib.L2P23_LINE
L2P3_LINE, # _xraylib.L2P3_LINE
L2P4_LINE, # _xraylib.L2P4_LINE
L2P5_LINE, # _xraylib.L2P5_LINE
L2Q1_LINE, # _xraylib.L2Q1_LINE
L3M1_LINE, # _xraylib.L3M1_LINE
L3M2_LINE, # _xraylib.L3M2_LINE
L3M3_LINE, # _xraylib.L3M3_LINE
L3M4_LINE, # _xraylib.L3M4_LINE
L3M5_LINE, # _xraylib.L3M5_LINE
L3N1_LINE, # _xraylib.L3N1_LINE
L3N2_LINE, # _xraylib.L3N2_LINE
L3N3_LINE, # _xraylib.L3N3_LINE
L3N4_LINE, # _xraylib.L3N4_LINE
L3N5_LINE, # _xraylib.L3N5_LINE
L3N6_LINE, # _xraylib.L3N6_LINE
L3N7_LINE, # _xraylib.L3N7_LINE
L3O1_LINE, # _xraylib.L3O1_LINE
L3O2_LINE, # _xraylib.L3O2_LINE
L3O3_LINE, # _xraylib.L3O3_LINE
L3O4_LINE, # _xraylib.L3O4_LINE
L3O45_LINE, # _xraylib.L3O45_LINE
L3O5_LINE, # _xraylib.L3O5_LINE
L3O6_LINE, # _xraylib.L3O6_LINE
L3O7_LINE, # _xraylib.L3O7_LINE
L3P1_LINE, # _xraylib.L3P1_LINE
L3P2_LINE, # _xraylib.L3P2_LINE
L3P23_LINE, # _xraylib.L3P23_LINE
L3P3_LINE, # _xraylib.L3P3_LINE
L3P4_LINE, # _xraylib.L3P4_LINE
L3P45_LINE, # _xraylib.L3P45_LINE
L3P5_LINE, # _xraylib.L3P5_LINE
L3Q1_LINE]

M_LINES = [
M1M2_LINE, # _xraylib.M1M2_LINE
M1M3_LINE, # _xraylib.M1M3_LINE
M1M4_LINE, # _xraylib.M1M4_LINE
M1M5_LINE, # _xraylib.M1M5_LINE
M1N1_LINE, # _xraylib.M1N1_LINE
M1N2_LINE, # _xraylib.M1N2_LINE
M1N3_LINE, # _xraylib.M1N3_LINE
M1N4_LINE, # _xraylib.M1N4_LINE
M1N5_LINE, # _xraylib.M1N5_LINE
M1N6_LINE, # _xraylib.M1N6_LINE
M1N7_LINE, # _xraylib.M1N7_LINE
M1O1_LINE, # _xraylib.M1O1_LINE
M1O2_LINE, # _xraylib.M1O2_LINE
M1O3_LINE, # _xraylib.M1O3_LINE
M1O4_LINE, # _xraylib.M1O4_LINE
M1O5_LINE, # _xraylib.M1O5_LINE
M1O6_LINE, # _xraylib.M1O6_LINE
M1O7_LINE, # _xraylib.M1O7_LINE
M1P1_LINE, # _xraylib.M1P1_LINE
M1P2_LINE, # _xraylib.M1P2_LINE
M1P3_LINE, # _xraylib.M1P3_LINE
M1P4_LINE, # _xraylib.M1P4_LINE
M1P5_LINE, # _xraylib.M1P5_LINE
M2M3_LINE, # _xraylib.M2M3_LINE
M2M4_LINE, # _xraylib.M2M4_LINE
M2M5_LINE, # _xraylib.M2M5_LINE
M2N1_LINE, # _xraylib.M2N1_LINE
M2N2_LINE, # _xraylib.M2N2_LINE
M2N3_LINE, # _xraylib.M2N3_LINE
M2N4_LINE, # _xraylib.M2N4_LINE
M2N5_LINE, # _xraylib.M2N5_LINE
M2N6_LINE, # _xraylib.M2N6_LINE
M2N7_LINE, # _xraylib.M2N7_LINE
M2O1_LINE, # _xraylib.M2O1_LINE
M2O2_LINE, # _xraylib.M2O2_LINE
M2O3_LINE, # _xraylib.M2O3_LINE
M2O4_LINE, # _xraylib.M2O4_LINE
M2O5_LINE, # _xraylib.M2O5_LINE
M2O6_LINE, # _xraylib.M2O6_LINE
M2O7_LINE, # _xraylib.M2O7_LINE
M2P1_LINE, # _xraylib.M2P1_LINE
M2P2_LINE, # _xraylib.M2P2_LINE
M2P3_LINE, # _xraylib.M2P3_LINE
M2P4_LINE, # _xraylib.M2P4_LINE
M2P5_LINE, # _xraylib.M2P5_LINE
M3M4_LINE, # _xraylib.M3M4_LINE
M3M5_LINE, # _xraylib.M3M5_LINE
M3N1_LINE, # _xraylib.M3N1_LINE
M3N2_LINE, # _xraylib.M3N2_LINE
M3N3_LINE, # _xraylib.M3N3_LINE
M3N4_LINE, # _xraylib.M3N4_LINE
M3N5_LINE, # _xraylib.M3N5_LINE
M3N6_LINE, # _xraylib.M3N6_LINE
M3N7_LINE, # _xraylib.M3N7_LINE
M3O1_LINE, # _xraylib.M3O1_LINE
M3O2_LINE, # _xraylib.M3O2_LINE
M3O3_LINE, # _xraylib.M3O3_LINE
M3O4_LINE, # _xraylib.M3O4_LINE
M3O5_LINE, # _xraylib.M3O5_LINE
M3O6_LINE, # _xraylib.M3O6_LINE
M3O7_LINE, # _xraylib.M3O7_LINE
M3P1_LINE, # _xraylib.M3P1_LINE
M3P2_LINE, # _xraylib.M3P2_LINE
M3P3_LINE, # _xraylib.M3P3_LINE
M3P4_LINE, # _xraylib.M3P4_LINE
M3P5_LINE, # _xraylib.M3P5_LINE
M3Q1_LINE, # _xraylib.M3Q1_LINE
M4M5_LINE, # _xraylib.M4M5_LINE
M4N1_LINE, # _xraylib.M4N1_LINE
M4N2_LINE, # _xraylib.M4N2_LINE
M4N3_LINE, # _xraylib.M4N3_LINE
M4N4_LINE, # _xraylib.M4N4_LINE
M4N5_LINE, # _xraylib.M4N5_LINE
M4N6_LINE, # _xraylib.M4N6_LINE
M4N7_LINE, # _xraylib.M4N7_LINE
M4O1_LINE, # _xraylib.M4O1_LINE
M4O2_LINE, # _xraylib.M4O2_LINE
M4O3_LINE, # _xraylib.M4O3_LINE
M4O4_LINE, # _xraylib.M4O4_LINE
M4O5_LINE, # _xraylib.M4O5_LINE
M4O6_LINE, # _xraylib.M4O6_LINE
M4O7_LINE, # _xraylib.M4O7_LINE
M4P1_LINE, # _xraylib.M4P1_LINE
M4P2_LINE, # _xraylib.M4P2_LINE
M4P3_LINE, # _xraylib.M4P3_LINE
M4P4_LINE, # _xraylib.M4P4_LINE
M4P5_LINE, # _xraylib.M4P5_LINE
M5N1_LINE, # _xraylib.M5N1_LINE
M5N2_LINE, # _xraylib.M5N2_LINE
M5N3_LINE, # _xraylib.M5N3_LINE
M5N4_LINE, # _xraylib.M5N4_LINE
M5N5_LINE, # _xraylib.M5N5_LINE
M5N6_LINE, # _xraylib.M5N6_LINE
M5N7_LINE, # _xraylib.M5N7_LINE
M5O1_LINE, # _xraylib.M5O1_LINE
M5O2_LINE, # _xraylib.M5O2_LINE
M5O3_LINE, # _xraylib.M5O3_LINE
M5O4_LINE, # _xraylib.M5O4_LINE
M5O5_LINE, # _xraylib.M5O5_LINE
M5O6_LINE, # _xraylib.M5O6_LINE
M5O7_LINE, # _xraylib.M5O7_LINE
M5P1_LINE, # _xraylib.M5P1_LINE
M5P2_LINE, # _xraylib.M5P2_LINE
M5P3_LINE, # _xraylib.M5P3_LINE
M5P4_LINE, # _xraylib.M5P4_LINE
M5P5_LINE]

N_LINES = [
N1N2_LINE, # _xraylib.N1N2_LINE
N1N3_LINE, # _xraylib.N1N3_LINE
N1N4_LINE, # _xraylib.N1N4_LINE
N1N5_LINE, # _xraylib.N1N5_LINE
N1N6_LINE, # _xraylib.N1N6_LINE
N1N7_LINE, # _xraylib.N1N7_LINE
N1O1_LINE, # _xraylib.N1O1_LINE
N1O2_LINE, # _xraylib.N1O2_LINE
N1O3_LINE, # _xraylib.N1O3_LINE
N1O4_LINE, # _xraylib.N1O4_LINE
N1O5_LINE, # _xraylib.N1O5_LINE
N1O6_LINE, # _xraylib.N1O6_LINE
N1O7_LINE, # _xraylib.N1O7_LINE
N1P1_LINE, # _xraylib.N1P1_LINE
N1P2_LINE, # _xraylib.N1P2_LINE
N1P3_LINE, # _xraylib.N1P3_LINE
N1P4_LINE, # _xraylib.N1P4_LINE
N1P5_LINE, # _xraylib.N1P5_LINE
N2N3_LINE, # _xraylib.N2N3_LINE
N2N4_LINE, # _xraylib.N2N4_LINE
N2N5_LINE, # _xraylib.N2N5_LINE
N2N6_LINE, # _xraylib.N2N6_LINE
N2N7_LINE, # _xraylib.N2N7_LINE
N2O1_LINE, # _xraylib.N2O1_LINE
N2O2_LINE, # _xraylib.N2O2_LINE
N2O3_LINE, # _xraylib.N2O3_LINE
N2O4_LINE, # _xraylib.N2O4_LINE
N2O5_LINE, # _xraylib.N2O5_LINE
N2O6_LINE, # _xraylib.N2O6_LINE
N2O7_LINE, # _xraylib.N2O7_LINE
N2P1_LINE, # _xraylib.N2P1_LINE
N2P2_LINE, # _xraylib.N2P2_LINE
N2P3_LINE, # _xraylib.N2P3_LINE
N2P4_LINE, # _xraylib.N2P4_LINE
N2P5_LINE, # _xraylib.N2P5_LINE
N3N4_LINE, # _xraylib.N3N4_LINE
N3N5_LINE, # _xraylib.N3N5_LINE
N3N6_LINE, # _xraylib.N3N6_LINE
N3N7_LINE, # _xraylib.N3N7_LINE
N3O1_LINE, # _xraylib.N3O1_LINE
N3O2_LINE, # _xraylib.N3O2_LINE
N3O3_LINE, # _xraylib.N3O3_LINE
N3O4_LINE, # _xraylib.N3O4_LINE
N3O5_LINE, # _xraylib.N3O5_LINE
N3O6_LINE, # _xraylib.N3O6_LINE
N3O7_LINE, # _xraylib.N3O7_LINE
N3P1_LINE, # _xraylib.N3P1_LINE
N3P2_LINE, # _xraylib.N3P2_LINE
N3P3_LINE, # _xraylib.N3P3_LINE
N3P4_LINE, # _xraylib.N3P4_LINE
N3P5_LINE, # _xraylib.N3P5_LINE
N4N5_LINE, # _xraylib.N4N5_LINE
N4N6_LINE, # _xraylib.N4N6_LINE
N4N7_LINE, # _xraylib.N4N7_LINE
N4O1_LINE, # _xraylib.N4O1_LINE
N4O2_LINE, # _xraylib.N4O2_LINE
N4O3_LINE, # _xraylib.N4O3_LINE
N4O4_LINE, # _xraylib.N4O4_LINE
N4O5_LINE, # _xraylib.N4O5_LINE
N4O6_LINE, # _xraylib.N4O6_LINE
N4O7_LINE, # _xraylib.N4O7_LINE
N4P1_LINE, # _xraylib.N4P1_LINE
N4P2_LINE, # _xraylib.N4P2_LINE
N4P3_LINE, # _xraylib.N4P3_LINE
N4P4_LINE, # _xraylib.N4P4_LINE
N4P5_LINE, # _xraylib.N4P5_LINE
N5N6_LINE, # _xraylib.N5N6_LINE
N5N7_LINE, # _xraylib.N5N7_LINE
N5O1_LINE, # _xraylib.N5O1_LINE
N5O2_LINE, # _xraylib.N5O2_LINE
N5O3_LINE, # _xraylib.N5O3_LINE
N5O4_LINE, # _xraylib.N5O4_LINE
N5O5_LINE, # _xraylib.N5O5_LINE
N5O6_LINE, # _xraylib.N5O6_LINE
N5O7_LINE, # _xraylib.N5O7_LINE
N5P1_LINE, # _xraylib.N5P1_LINE
N5P2_LINE, # _xraylib.N5P2_LINE
N5P3_LINE, # _xraylib.N5P3_LINE
N5P4_LINE, # _xraylib.N5P4_LINE
N5P5_LINE, # _xraylib.N5P5_LINE
N6N7_LINE, # _xraylib.N6N7_LINE
N6O1_LINE, # _xraylib.N6O1_LINE
N6O2_LINE, # _xraylib.N6O2_LINE
N6O3_LINE, # _xraylib.N6O3_LINE
N6O4_LINE, # _xraylib.N6O4_LINE
N6O5_LINE, # _xraylib.N6O5_LINE
N6O6_LINE, # _xraylib.N6O6_LINE
N6O7_LINE, # _xraylib.N6O7_LINE
N6P1_LINE, # _xraylib.N6P1_LINE
N6P2_LINE, # _xraylib.N6P2_LINE
N6P3_LINE, # _xraylib.N6P3_LINE
N6P4_LINE, # _xraylib.N6P4_LINE
N6P5_LINE, # _xraylib.N6P5_LINE
N7O1_LINE, # _xraylib.N7O1_LINE
N7O2_LINE, # _xraylib.N7O2_LINE
N7O3_LINE, # _xraylib.N7O3_LINE
N7O4_LINE, # _xraylib.N7O4_LINE
N7O5_LINE, # _xraylib.N7O5_LINE
N7O6_LINE, # _xraylib.N7O6_LINE
N7O7_LINE, # _xraylib.N7O7_LINE
N7P1_LINE, # _xraylib.N7P1_LINE
N7P2_LINE, # _xraylib.N7P2_LINE
N7P3_LINE, # _xraylib.N7P3_LINE
N7P4_LINE, # _xraylib.N7P4_LINE
N7P5_LINE]