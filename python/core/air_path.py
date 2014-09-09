from sim_class import *
O2 = compound(CP='O2')
N2 = compound(CP='N2')
CO2 = compound(CP='CO2')
H2 = compound(CP='H2')
Ar = compound(CP='Ar')
Ne = compound(CP='Ne')
He = compound(CP='He')
Kr = compound(CP='Kr')
Xe = compound(CP='Xe')

AIR = mix_compound([O2, N2, CO2, H2, Ar, Ne, He, Kr, Xe], 
	[23.2,75.47,.0046,0,1.28,.0012,.00007,.0003,.00004])
	
def air_atten(ev, y, l, gas = AIR, density = 0.001225):
	return y*np.exp(-gas.mac_total(ev)*density*l)