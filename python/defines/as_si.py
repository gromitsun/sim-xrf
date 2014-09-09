from sim_class import *

RECALC = False
calc_range = [0,180,0,90]
# calc_range = [180-5, 180, -90, 90]
omega_calc = solid_angle(angle_range = calc_range)
datapath = 'D:/desktop/data/'
filename = 'as_si_500um_psi45.h5'

# sample definition
As = compound(CP = 'As')
As = sample(As, element_density(symbol2number('As')),5e-8)

Si = compound(CP = 'Si')
Si = sample(Si,element_density(symbol2number('Si')),500e-4) ############

s = sample_multilayer([As,Si])

# detector definition
det_center = 90
x = 10
y = 10
if det_center == 90:
	angle_range_gen = lambda x: [90-x, 90+x, -x, x]
elif det_center == 180:
	angle_range_gen = lambda x: [180-x, 180, -90, 90]
elif det_center == 0:
	angle_range_gen = lambda x: [0, x, 0, 360]
_omega=solid_angle(angle_range = angle_range_gen(x))
ch = channel()
resp = response()
det = detector(ch, resp, omega = deepcopy(_omega))

# illumination definition
il = illumination(ev0=1.2e4, psi=45, alpha=0)

# plot specs
PLOT_SPEC = True
DET_RESPONSE = True
xlim = [1e2,1.3e4]
ylim = [1e-18,5e-5]
# peak to background and snr calculations
CALC_SNR = True
from sim_snr import *
peak_center = 10540 # ev
bg_leftshift = 500 # ev
bg_rightshift = 500 # ev

pnb_func = pnb_fixedwidth
pnb_func = pnb_ascalc
pnb_func = pnb_expfit
pnb_func = [pnb_fixedwidth, pnb_expfit, pnb_ascalc]
pnb_func = [pnb_fixedwidth, pnb_expfit, pnb_ascalc, pnb_snip]
element = 'As' # needed for pnb_ascalc
openings = np.arange(1,91)