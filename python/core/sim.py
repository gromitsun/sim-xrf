from sim_calcspec import *
from sim_genspec import *
import os
import sys
from sim_postprocess import *
from sim_snr import *
sys.path.append(os.path.realpath('../defines/'))

argv_list = sys.argv[1:]


def beta_symmetry(beta_range, axis0 = True, axis90 = True):
	if axis0:
		beta_range = [abs(x) for x in beta_range]
	beta_range = [(x+180)%360-180 for x in beta_range]
	if axis90 and abs(beta_range[0])>=90 and abs(beta_range[1])>=90:
		beta_range = [180-x for x in beta_range]
	beta_range = [(x+180)%360-180 for x in beta_range]
	beta_range.sort()
	return beta_range

def main(def_file, show = True, block = True):
	try:
		exec 'import '+def_file+' as define'
	except ImportError:
		print 'Cannot import \''+def_file+'\'.'
		import sim_define as define
	
	RECALC = define.RECALC
	datapath = define.datapath
	filename = define.filename
	s = define.s
	omega_calc = define.omega_calc
	det = define.det
	il = define.il
	DET_RESPONSE = define.DET_RESPONSE
	PLOT_SPEC = define.PLOT_SPEC
	CALC_SNR = define.CALC_SNR
	det_center = define.det_center
	CALC_SNR = define.CALC_SNR
	xlim = define.xlim
	ylim = define.ylim
	openings = define.openings
	pnb_func = define.pnb_func
	angle_range_gen = define.angle_range_gen
	element = define.element
	peak_center = define.peak_center
	bg_leftshift = define.bg_leftshift
	bg_rightshift = define.bg_rightshift
	
	try:
		snr_kwargs = define.snr_kwargs
	except AttributeError:
		snr_kwargs = {}
		
	
	if RECALC or not os.path.exists(datapath+filename):
		# calculate spectrum
		if s:
			xrf, ray, comp = rawspec_multilayer(s,il,omega_calc)
			ray = [sum_spec(ray,mat=True)]
			comp = [sum_spec(comp,mat=True)]
		else:
			compound = define.c
			xrf = rawspec_xrf(compound, il.ev0, omega=omega_calc)
			ray = [rawspec_rayleigh(compound, il.ev0, omega=omega_calc)]
			comp = [rawspec_compton(compound, il.ev0, omega=omega_calc)]

		# save spectrum
		savespec(xrf+ray+comp, filename = filename, dir = datapath)	

	# read spectrum from file
	rawspec_arr = readspec(filename = filename, dir = datapath)	

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
	if 	PLOT_SPEC:
		plt.figure(filename[:-3]+'_theta'+str(det_center)+'_spec')
		if not DET_RESPONSE:
			yshift = ylim[0]
		else:
			yshift = 0
		plot_specs(rawspec_arr,det, ylim = ylim, xlim = xlim, show = False, show_total = DET_RESPONSE, yshift = yshift,save_npz=True, npz_name='spec_'+filename[:-3]+'_theta'+str(det_center)+'_sa15.npz')

		
	# SNR calculation
	if CALC_SNR:
		total_arr = []
		signal_arr = []
		for x in openings:
			omega_slice = omega_calc.slice(angle_range = angle_range_gen(x))
			for rawspec in rawspec_arr:
				ev_arr, intensity_arr = genspec_raw(rawspec.ev_mat, rawspec.intensity_mat*omega_slice)
				rawspec.ev_arr, rawspec.intensity_arr = genspec_det(ev_arr, intensity_arr, det)
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
			ax1.plot(openings,SNR, 'r-', label = 'S/N')
			# from matplotlib.ticker import ScalarFormatter as sfmt
			# sfmt.set_powerlimits((-3,3))
			# sfmt.set_scientific(True)
			ax1.ticklabel_format(style='sci', axis='y', scilimits=(-2,2))
			plt.ylabel('S/N (red solid)', color='r')
			plt.xlabel('opening angle (deg)')
			for tl in ax1.get_yticklabels():
				tl.set_color('r')
			ax2 = ax1.twinx()
			ax2.plot(openings,p2b, 'b--', label = 'P/B')
			ax2.ticklabel_format(style='sci', axis='y', scilimits=(-2,2))
			plt.ylabel('P/B (blue dashed)',color = 'b')
			for tl in ax2.get_yticklabels():
				tl.set_color('b')
			
			# # second x axis
			# xmin, xmax = ax1.xaxis.get_view_interval()
			# ax3 = ax2.twiny()
			# # Make some room at the bottom
			# fig.subplots_adjust(bottom=0.2)
			# ax3.set_frame_on(True)
			# ax3.patch.set_visible(False)
			# ax3.xaxis.set_ticks_position('bottom')
			# ax3.xaxis.set_label_position('bottom')
			# ax3.spines['bottom'].set_position(('outward', 50))
			# xmin, xmax = [solid_angle(angle_range = angle_range_gen(x)).subtend for x in [xmin,xmax]]
			# ax3.set_xlim(xmin=xmin, xmax=xmax)
			# tick_locs = np.arange(xmin,int(xmax)+1)
			# plt.xticks(tick_locs)
			# plt.xlabel('subtended angle (sr)')
		np.savez('snr_wpw'+str(det_center)+filename[19:-9]+'.npz',snr=SNR,p2b=p2b,openings=openings)
		# np.savez('snr_as'+str(det_center)+filename[6:-9]+'.npz',snr=SNR,p2b=p2b,openings=openings)
	# show all plots
	if show:
		plt.show(block = block)

for def_file in argv_list:
	if def_file == argv_list[-1]:
		block = True
	else:
		block = False
	main(def_file, show = True, block = block)
	


# # # # # # multiprocessing.Pool
	
# def run():
	# for def_file in argv_list:
		# main(def_file)

# from multiprocessing import Pool, freeze_support

# if __name__ == '__main__':
	# freeze_support()
	# pool = Pool()
	# pool.apply_async(run)
	# pool.close()
	# pool.join()
# # # # # # # # # # # # # # # # # # # # #


# # # # # # # # # # # # multiprocessing.Process

# from multiprocessing import Process, freeze_support

# if __name__ == '__main__':
	# freeze_support()
	# for def_file in argv_list:
		# p = Process(target = main, args = (def_file,))
		# p.start()
		# p.join()