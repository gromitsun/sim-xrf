from sim_class import *

RECALC = False
calc_range = [0,180,0,90]
# calc_range = [180-5, 180, -90, 90]
omega_calc = solid_angle(angle_range = calc_range)
datapath = 'D:/desktop/data/'
filename = 'water_protein_water_100nm_psi90.h5'

# sample definition
protein = compound(CP = 'C30H50O10N9S')
water = compound(CP = 'H2O')
Zn = compound(CP = 'Zn')

p1 = mix_compound([protein, Zn],[1,1e-4])
protein = sample(p1,1.05,20e-7)

water = sample(water,1,100e-7)

s = sample_multilayer([water, protein, deepcopy(water)])

# detector definition
det_center = 180
x = 15
y = 60
if det_center == 90:
	angle_range_gen = lambda x: [90-x, 90+x, -x, x]
elif det_center == 180:
	angle_range_gen = lambda x: [180-x, 180, -90, 90]
elif det_center == 0:
	angle_range_gen = lambda x: [0, x, 0, 360]
_omega=solid_angle(angle_range = angle_range_gen(x))
ch = channel(0.,10.,2048)
resp = response()
det = detector(ch, resp, omega = deepcopy(_omega))

# illumination definition
il = illumination(ev0=1e4, psi=90, alpha=0)

# plot specs
PLOT_SPEC = True
DET_RESPONSE = True
xlim = [1e2,1.1e4]
xlim = [0,11]
ylim = [1e-18,5e-6]
ylim = [1e-18,7e-8]

# peak to background and snr calculations
CALC_SNR = False
from sim_snr import *
peak_center = 8630 # ev
bg_leftshift = 500 # ev
bg_rightshift = 500 # ev
pnb_func = pnb_expfit
pnb_func = pnb_fixedwidth
pnb_func = [pnb_fixedwidth, pnb_expfit, pnb_ascalc, pnb_snip]
pnb_func = pnb_snip
pnb_func = pnb_ascalc
snr_kwargs = {'loops':8, 'end_loops':6}
element = 'Zn' # needed for pnb_ascalc
openings = np.arange(1,91)