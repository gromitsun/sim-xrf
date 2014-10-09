import numpy as np

def pnb_fixedwidth(ev, y, peak_center, bg_range, peak_width = 100):
	
	def ev2ch(x):
		return np.abs(ev-x).argmin()
	
	b = np.array([ev2ch(x) for x in bg_range])
	p = np.array([ev2ch(peak_center-peak_width/2.),ev2ch(peak_center+peak_width/2.)])
	
	Nt = np.sum(y[p[0]:p[1]])
	B = 0
	nchB = 0
	for i in range(0,len(b),2):
		B += np.sum(y[b[i]:b[i+1]])
		nchB = b[i+1] - b[i]
	B *= 1.*(p[1] - p[0])/nchB
	return Nt-B, B