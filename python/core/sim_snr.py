import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

from sim_common import flatten

import matplotlib.pyplot as plt

def listgen(spec_list,select = None):
	"""
	list of 'spectrum' objects -> 2D numpy array of stacked spectra with name = select
	"""
	if select != None:
		_spec_list = []
		for spec in flatten(spec_list):
			if spec.name == select:
				_spec_list.append(spec)
		spec_list = _spec_list
	return np.array([spec.intensity_arr for spec in spec_list])
			

def pnb_ascalc(ev, total, spec, peak_center, bg_leftshift, bg_rightshift, detector, peak_width = 1.0, **kwargs):
	b1 = peak_center - bg_leftshift # ev
	b2 = peak_center + bg_rightshift # ev
	
	fano = detector.response.fano
	noise = detector.response.noise # FWHM
	ev2ch = detector.channel.channel
	FWHM =  detector.response.FWHM(peak_center)*peak_width
	
	nb1 = int(ev2ch(b1))
	nb2 = int(ev2ch(b2))

	npeak1 = int(ev2ch(peak_center-FWHM))
	npeak2 = int(ev2ch(peak_center+FWHM))
	
	Nt = np.sum(total[:,npeak1:npeak2],axis=1)
	P = np.sum(spec[:,npeak1:npeak2],axis=1)
	return P, Nt-P

def pnb_fixedwidth(ev, total, peak_center, bg_leftshift, bg_rightshift, detector, peak_width = 1.0, **kwargs):
	bg_leftshift = 250
	bg_rightshift = 250
	
	b1 = peak_center - bg_leftshift # ev
	b2 = peak_center + bg_rightshift # ev

	fano = detector.response.fano
	noise = detector.response.noise # FWHM
	ch_width = detector.channel.gain
	ev2ch = detector.channel.channel
	sigma = np.sqrt((noise/2.3548)**2+3.58*fano*peak_center)
	FWHM = 2.3548*sigma
	
	
	
	nb1 = int(ev2ch(b1))
	nb2 = int(ev2ch(b2))
	ncenter = int(ev2ch(peak_center))
	
	nhalfwidth = int(np.round(FWHM/ch_width))
	
	nb1 -= nhalfwidth/2
	nb2 -= nhalfwidth/2
	
	Nt = np.sum(total[:,ncenter-nhalfwidth:ncenter+nhalfwidth],axis=1)
	B = np.sum(total[:,nb1:nb1+nhalfwidth],axis=1)+np.sum(total[:,nb2:nb2+nhalfwidth],axis=1)
	return Nt-B, B
	
	
def pnb_interp(ev, total, peak_center, bg_leftshift, bg_rightshift, **kwargs):
	b1 = peak_center - bg_leftshift # ev
	b2 = peak_center + bg_rightshift # ev
	
	nb1 = np.where(ev == b1)[0][0]
	nb2 = np.where(ev == b2)[0][0]
	ncenter = np.where(ev == peak_center)[0][0]

	Nt = np.array([])
	B = np.array([])

	for i in range(len(total)):
		x = np.append(ev[nb1-3:nb1+3],ev[nb2-3:nb2+3])
		y = np.append(total[i,nb1-3:nb1+3],total[i,nb2-3:nb2+3])
		background = np.exp(interp1d(x,np.log(y),kind='linear')(ev[nb1:nb2+1]))
		signal = total[i,nb1:nb2+1]
		nrange = np.where(signal>background)
		
		
		
		# # # show interp results
		if i>9:
			plt.plot(ev[nb1:nb2+1],signal)
			plt.plot(ev[nb1:nb2+1],background)
			plt.show()
			raw_input()
		
		if len(nrange[0]):
			Nt = np.append(Nt,signal[nrange].sum())
			B = np.append(B,background[nrange].sum())
		else:
			Nt = np.append(Nt,1)
			B = np.append(B,1)
	return Nt-B, B
	
	
def pnb_expfit(ev, total, peak_center, bg_leftshift, bg_rightshift, detector, peak_width = 1.0, **kwargs):
	bg_leftshift = 250
	bg_rightshift = 250
	
	b1 = peak_center - bg_leftshift # ev
	b2 = peak_center + bg_rightshift # ev
	
	fano = detector.response.fano
	noise = detector.response.noise # FWHM
	ch_width = detector.channel.gain
	ev2ch = detector.channel.channel
	FWHM = detector.response.FWHM(peak_center)*peak_width
	
	nb1 = int(ev2ch(b1))
	nb2 = int(ev2ch(b2))
	ncenter = int(ev2ch(peak_center))
	
	nhalfwidth = int(np.round(FWHM/ch_width))
	
	nb10 = nb1-nhalfwidth/2
	nb20 = nb2-nhalfwidth/2
	
	Nt = np.array([])
	B = np.array([])


	def exponential(x,b,c,d):
		return np.exp(b*(x+c))+d
	def polynomial(x,a,b,c,d):
		return a+b*x+c*x**2+d*x**3
	func = polynomial
	
	for i in range(len(total)):
		x = np.append(ev[nb10:nb10+nhalfwidth],ev[nb20:nb20+nhalfwidth])
		y = np.append(total[i,nb10:nb10+nhalfwidth],total[i,nb20:nb20+nhalfwidth])

	
		# curve fit
		x1, x2 = (ev[nb1], ev[nb2])
		y1, y2 = (total[i,nb1], total[i,nb2])
		a0 = 1.*(x1*y2-x2*y1)/(x1-x2)
		b0 = 1.*(y2-y1)/(x2-x1)
		c0 = 0
		d0 = 0
		popt, pcov = curve_fit(func, x, y, p0=[a0,b0,c0,d0])
		# print popt
		
		background = func(ev[ncenter-nhalfwidth:ncenter+nhalfwidth],*popt)
		signal = total[i,ncenter-nhalfwidth:ncenter+nhalfwidth]
		# # # show interp results
		# """
		if i in [4,44,89]:
			plt.figure(str(i+1))
			plt.plot(ev[nb10:nb20+nhalfwidth],total[i,nb10:nb20+nhalfwidth])
			plt.plot(ev[nb10:nb20+nhalfwidth],func(ev[nb10:nb20+nhalfwidth],*popt))
			# plt.show()
			# raw_input()
		# """
			
		Nt = np.append(Nt,signal.sum())
		B = np.append(B,background.sum())
	return Nt-B, B
	
import csnip as snip
def pnb_snip(ev, total, peak_center, bg_leftshift, bg_rightshift, detector, peak_width = 1.0, **kwargs):
	b1 = peak_center - bg_leftshift # ev
	b2 = peak_center + bg_rightshift # ev

	fano = detector.response.fano
	noise = detector.response.noise # FWHM
	ch_width = detector.channel.gain
	ev2ch = detector.channel.channel
	FWHM_func =  detector.response.FWHM
	FWHM =  FWHM_func(peak_center)*peak_width
	
	nb1 = int(ev2ch(b1))
	nb2 = int(ev2ch(b2))

	npeak1 = int(ev2ch(peak_center-FWHM))
	npeak2 = int(ev2ch(peak_center+FWHM))
	
	# ev = ev[nb1:nb2+1].reshape(-1,1)
	# total = total[:,nb1:nb2+1].T
	
	# B = snip.snip(ev,total, FWHM = FWHM_func)
	
	# plt.figure()
	# plt.plot(ev,B[:,(4,44,89,)],'--')
	# plt.plot(ev,total[:,(4,44,89,)],'-')
	
	# Nt = total[npeak1-nb1:npeak2-nb1+1].sum(axis = 0)
	# B = B[npeak1-nb1:npeak2-nb1+1].sum(axis = 0)
	
	ev = ev[nb1:nb2+1]
	total = total[:,nb1:nb2+1]
	
	B = np.array([snip.snip(ev,y, detector = detector, **kwargs) for y in total])
	
	# # # show background clipping results
	for i in [4,44,89]:
		plt.figure(str(i+1))
		plt.plot(ev,total[i,:].T,'-')
		plt.plot(ev,B[i,:].T,'--')
	
	Nt = total[:,npeak1-nb1:npeak2-nb1+1].sum(axis = 1)
	B = B[:,npeak1-nb1:npeak2-nb1+1].sum(axis = 1)

	return Nt-B, B