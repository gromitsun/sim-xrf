from sim_class import *

RECALC = False
calc_range = [0,180,0,90]
# calc_range = [180-5, 180, -90, 90]
omega_calc = solid_angle(angle_range = calc_range)
datapath = 'D:/desktop/data/'
filename = 'ti_si_1nm_1um_psi45.h5'

# sample definition
Ti = compound(CP = 'Ti')
Ti = sample(Ti, element_density(symbol2number('Ti')),1e-8)

Si = compound(CP = 'Si')
Si = sample(Si, element_density(symbol2number('Si')),1e-1)

Si3N4 = compound(CP = 'Si3N4')
Si3N4 = sample(Si3N4,3.44,1e-4) ############

s = sample_multilayer([Ti,Si])

# detector definition
det_center = 90
x = 5
y = 5
if det_center == 90:
	angle_range_gen = lambda x: [90-x, 90+x, -x, x]
elif det_center == 180:
	angle_range_gen = lambda x: [180-x, 180, -90, 90]
elif det_center == 0:
	angle_range_gen = lambda x: [0, x, 0, 360]
_omega=solid_angle(angle_range = angle_range_gen(x))
ch = channel(1e2,1.2e4,10)
resp = response()
det = detector(ch, resp, omega = deepcopy(_omega))

# illumination definition
il = illumination(ev0=1e4, psi=45, alpha=0)

# plot specs
PLOT_SPEC = True
xlim = [1e2,1.3e4]
ylim = [1e-18,5e-6]
# peak to background and snr calculations
CALC_SNR = True
from sim_snr import *
peak_center = 4510 # ev
bg_leftshift = 250 # ev
bg_rightshift = 250 # ev

pnb_func = pnb_fixedwidth
pnb_func = pnb_ascalc
pnb_func = pnb_expfit
element = 'Ti' # needed for pnb_ascalc
openings = np.arange(1,91)