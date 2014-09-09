import sys
sys.path.append('../../core/')

from sim_calcspec import *
from sim_genspec import *
from sim_postprocess import *
from sim_snr import *

filename = 'as_t_n.h5'
x = 6.7
calc_range = [90-x, 90+x, -x, x]
# calc_range = [180-52.6, 180-13.9, 0, 90]
omega_calc = solid_angle(angle_range = calc_range)

# sample definition
As = compound(CP = 'As')
s_As = sample(As, element_density(symbol2number('As')),5e-8)

Si = compound(CP = 'Si')
s_Si = sample(Si,element_density(symbol2number('Si')),500e-4) ############

s = sample_multilayer([s_As,s_Si])

# detector definition
det_center = 90
x = 6.7
y = 6.7
if det_center == 90:
	angle_range_gen = lambda x: [90-x, 90+x, -x, x]
elif det_center == 180:
	angle_range_gen = lambda x: [180-52.6, 180-13.9, 0, 360]
elif det_center == 0:
	angle_range_gen = lambda x: [0, x, 0, 360]
_omega=solid_angle(angle_range = angle_range_gen(x))
ch = channel(0.,10.,2048)
resp = response()
det = detector(ch, resp, omega = deepcopy(_omega))

# illumination definition
il = illumination(ev0=1.2e4, psi=45, alpha=0)

# plot specs
xlim = [0,11]
ylim = [1e-18,5e-6]
DET_RESPONSE = True

# peak to background and snr calculations
peak_center = 10540 # ev
bg_leftshift = 500 # ev
bg_rightshift = 500 # ev
pnb_func = pnb_expfit
pnb_func = pnb_fixedwidth
pnb_func = pnb_snip
pnb_func = pnb_ascalc
snr_kwargs = {'loops':8, 'end_loops':6}
element = 'As' # needed for pnb_ascalc


def beta_symmetry(beta_range, axis0 = True, axis90 = True):
	if axis0:
		beta_range = [abs(x) for x in beta_range]
	beta_range = [(x+180)%360-180 for x in beta_range]
	if axis90 and abs(beta_range[0])>=90 and abs(beta_range[1])>=90:
		beta_range = [180-x for x in beta_range]
	beta_range = [(x+180)%360-180 for x in beta_range]
	beta_range.sort()
	return beta_range

total_arr = []
signal_arr = []
t_arr = np.linspace(500e-7,500e-4,100)
for t in t_arr:
	s_Si = sample(Si,element_density(symbol2number('Si')),t) ############
	s = sample_multilayer([s_As,s_Si])
	# calculate spectrum
	xrf, ray, comp = rawspec_multilayer(s,il,omega_calc)
	ray = [sum_spec(ray,mat=True)]
	comp = [sum_spec(comp,mat=True)]

	# read spectrum from file
	rawspec_arr = xrf+ray+comp

	# slice raw spectrum in the designated solid angle
	omega_calc = rawspec_arr[0].omega
	omega_slice = omega_calc.slice(angle_range=det.omega.range)


	if det.omega.range[0]<omega_calc.range[0] or det.omega.range[1]>omega_calc.range[1]:
		print 'Warning: theta out of range!'

	axis90 = (il.psi == 90)
	axis0 = (il.alpha == 0)

	if det.omega.range[2]<omega_calc.range[2]:
		beta_range = beta_symmetry([det.omega.range[2],omega_calc.range[2]-rad2deg(omega_calc.delta_beta)],axis0,axis90)
		omega_slice += omega_calc.slice(angle_range=det.omega.range[:2]+beta_range)
	elif det.omega.range[3]>omega_calc.range[3]:
		beta_range = beta_symmetry([omega_calc.range[3]+rad2deg(omega_calc.delta_beta),det.omega.range[3]],axis0,axis90)
		omega_slice += omega_calc.slice(angle_range=det.omega.range[:2]+beta_range)

	# generate spectrum w/ detector response
	for x in rawspec_arr:
		ev_arr, intensity_arr = genspec_raw(x.ev_mat, x.intensity_mat*omega_slice)
		if DET_RESPONSE:
			x.ev_arr, x.intensity_arr = genspec_det(ev_arr, intensity_arr, det)
		else:
			x.ev_arr, x.intensity_arr = genspec_bin(ev_arr, intensity_arr, det.channel.range)
		
	# plot spectrum
	if 	False:
		plt.figure(filename[:-3]+'_theta'+str(det_center)+'_spec')
		if not DET_RESPONSE:
			yshift = ylim[0]
		else:
			yshift = 0
		plot_specs(rawspec_arr,det, ylim = ylim, xlim = xlim, show = False, show_total = DET_RESPONSE, yshift = yshift,save_npz=True)

		
	# SNR calculation
	if True:
		_rawspec_total = sum_spec(rawspec_arr,arr=True)
		total_arr.append(_rawspec_total.intensity_arr)
		signal_arr.append(listgen(rawspec_arr,select = element).flatten())

		
total_arr = np.array(total_arr)
signal_arr = np.array(signal_arr)
ev_arr = _rawspec_total.ev_arr

for _pnb_func in flatten([pnb_func]):
	if _pnb_func == pnb_ascalc:
		P, B = _pnb_func(ev_arr, total_arr, signal_arr, peak_center, bg_leftshift, bg_rightshift, det, peak_width =1.)
	else:
		P, B = _pnb_func(ev_arr, total_arr, peak_center, bg_leftshift, bg_rightshift, det, peak_width =1.,**snr_kwargs)
	p2b = P/B
	SNR = P/np.sqrt(P+2*B)
	fig = plt.figure(filename[:-3]+'_theta'+str(det_center)+'_snr_'+_pnb_func.__name__)
	ax1 = fig.add_subplot(111)
	ax1.plot(t_arr*1e4,SNR, 'r-', label = 'S/N')
	ax1.ticklabel_format(style='sci', axis='y', scilimits=(-2,2))
	plt.ylabel('S/N (red solid)', color='r')
	plt.xlabel('water layer thickness (um)')
	for tl in ax1.get_yticklabels():
		tl.set_color('r')
	ax2 = ax1.twinx()
	ax2.plot(t_arr*1e4,p2b, 'b--', label = 'P/B')
	ax2.ticklabel_format(style='sci', axis='y', scilimits=(-2,2))
	plt.ylabel('P/B (blue dashed)',color = 'b')
	for tl in ax2.get_yticklabels():
		tl.set_color('b')
	plt.figure('Minimum required incident photons')
	plt.plot(t_arr*1e4,3./SNR)
	plt.ylabel(r'Minimum required incident photons ($\bar{n}$)')
	plt.xlabel('water layer thickness (um)')
	ax = plt.gca()
	ax.ticklabel_format(style='sci', axis='y', scilimits=(-2,2))

np.savez('as90.npz',snr = SNR, p2b = p2b, t_arr = t_arr)
# show all plots
plt.show()
