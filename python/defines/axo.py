from sim_class import *

RECALC = False
calc_range = [0,180,0,90]
# calc_range = [180-5, 180, -90, 90]
omega_calc = solid_angle(angle_range = calc_range)
datapath = 'D:/desktop/data/'
filename = 'axo_3nm_100nm_psi75.h5'

# sample definition
Pb = compound(CP = 'Pb')
La = compound(CP = 'La')
Pd = compound(CP = 'Pd')
Mo = compound(CP = 'Mo')
Cu = compound(CP = 'Cu')
Fe = compound(CP = 'Fe')
Ca = compound(CP = 'Ca')
SiN = compound(CP = 'Si3N4')

axo = mix_compound([Pb,La,Pd,Mo,Cu,Fe,Ca],
					[7.61,11.01,1.8,1.32,2.84,5.04,19.31])

s_axo = sample(axo,48.93/3,3e-7)
s_sin = sample(SiN,3.44,100e-7)

s = sample_multilayer([s_axo,s_sin])


# detector definition
det_center = 95
x = 6.7
y = 6.7
if det_center == 90:
	angle_range_gen = lambda x: [90-x, 90+x, -x, x]
elif det_center == 124:
	angle_range_gen = lambda x: [124-x, 124+x, -8, 8]
elif det_center == 50:
	angle_range_gen = lambda x: [50-x, 50+x, -8.9, 8.9]
elif det_center == 88:
	angle_range_gen = lambda x: [88-x, 88+x, -x, x]
elif det_center == 92:
	angle_range_gen = lambda x: [92-x, 92+x, -x, x]
elif det_center == 85:
	angle_range_gen = lambda x: [85-x, 85+x, -x, x]
elif det_center == 95:
	angle_range_gen = lambda x: [95-x, 95+x, -x, x]
_omega=solid_angle(angle_range = angle_range_gen(x))
ch = channel(2.97,12.26952,2048)
resp = response()
det = detector(ch, resp, omega = deepcopy(_omega))

# illumination definition
il = illumination(ev0=10.1e3, psi=75, alpha=0)

# plot specs
PLOT_SPEC = True
DET_RESPONSE = True
xlim = [0,11]
ylim = [1e-18,5e-6]

# peak to background and snr calculations
CALC_SNR = False
from sim_snr import *
peak_center = 8630 # ev
bg_leftshift = 500 # ev
bg_rightshift = 500 # ev
pnb_func = pnb_expfit
pnb_func = pnb_fixedwidth
pnb_func = pnb_ascalc
pnb_func = pnb_snip
pnb_func = [pnb_fixedwidth, pnb_expfit, pnb_ascalc, pnb_snip]
snr_kwargs = {'loops':8, 'end_loops':6}
element = 'Zn' # needed for pnb_ascalc
openings = np.arange(1,91)