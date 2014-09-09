from sim_class import *
from scipy.integrate import dblquad
x = 90
for x in range(1,91):
	omega=solid_angle(angle_range = [90-x, 90+x, -x, x], n_theta = 4*x, n_beta = 4*x)		
	x = np.radians(x)
	a,b,gfun,hfun = [np.pi/2-x, np.pi/2+x, lambda y: -x, lambda y: x]

	subtend = dblquad(lambda y,x: np.sin(x),a,b,gfun,hfun)[0]
	# subtend = dblquad(lambda y,x: x,0,1,lambda x:-1,lambda x:1)[0]

	print subtend,omega.subtend_calc,omega.subtend
