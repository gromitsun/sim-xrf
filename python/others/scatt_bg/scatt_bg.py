from sim_calcspec import *
import matplotlib.pyplot as plt
from scipy.integrate import dblquad
from matplotlib.colors import LogNorm

singleZ = 0
multiZ = 1
direct = 0
slice = 0
integrate = 1

ev = 1e4
Z = symbol2number('O')



if singleZ:
	if slice:
		x = 90
		y = 90
		omega=solid_angle(angle_range = [90-x, 90+x, -y, y], n_theta = 4*x, n_beta = 4*y)

		ray = dcs_mat_integral(dcs_rayleigh_pol,omega,ev,Z)
		comp = dcs_mat_integral(dcs_compton_pol,omega,ev,Z)
		ray_sum = np.array([])
		comp_sum = np.array([])

		for x in range(1,91):
			angle_range = [90-x, 90+x, -x, x]
			subtend = solid_angle(angle_range = angle_range, n_theta = 4*x, n_beta = 4*y).subtend
			theta_slice = np.logical_and(omega.theta_arr>=angle_range[0],omega.theta_arr<=angle_range[1])
			beta_slice = np.logical_and(omega.beta_arr>=angle_range[2],omega.beta_arr<=angle_range[3])
			omega_slice = np.multiply(*np.meshgrid(theta_slice,beta_slice))
			ray_sum = np.append(ray_sum,(omega_slice*ray).sum()/subtend**2)
			comp_sum = np.append(comp_sum,(omega_slice*comp).sum()/subtend**2)

	if direct:
		ray_sum = np.array([])
		comp_sum = np.array([])
		for x in range(1,91):
			omega=solid_angle(angle_range = [90-x, 90+x, -x, x], n_theta = 4*x, n_beta = 4*x)
			ray = dcs_mat_integral(dcs_rayleigh_pol,omega,ev,Z)/omega.subtend**2
			comp = dcs_mat_integral(dcs_compton_pol,omega,ev,Z)/omega.subtend**2
			ray_sum = np.append(ray_sum,ray.sum())
			comp_sum = np.append(comp_sum,comp.sum())
	if integrate:
		ray_sum = np.array([])
		comp_sum = np.array([])
		_dcs_rayleigh_pol = lambda beta_rad, theta_rad:dcs_rayleigh_pol(ev, Z, theta_rad, beta_rad)*np.sin(theta_rad)
		_dcs_compton_pol = lambda beta_rad, theta_rad:dcs_compton_pol(ev, Z, theta_rad, beta_rad)*np.sin(theta_rad)
		for x in range(1,91):
			omega=solid_angle(angle_range = [90-x, 90+x, -x, x], n_theta = 4*x, n_beta = 4*x)
			x = np.radians(x)
			a,b,gfun,hfun = [np.pi/2-x, np.pi/2+x, lambda y: -x, lambda y: x]
			ray = dblquad(_dcs_rayleigh_pol,a,b,gfun,hfun)[0]/omega.subtend**2
			comp = dblquad(_dcs_compton_pol,a,b,gfun,hfun)[0]/omega.subtend**2
			ray_sum = np.append(ray_sum,ray)
			comp_sum = np.append(comp_sum,comp)

	plt.figure(number2symbol(Z)+str(Z)+'_theta90')
	plt.title(number2symbol(Z)+' Z = %d, det@90$^\circ$'%Z)
	plt.plot(ray_sum, label = 'Rayleigh')
	plt.plot(comp_sum, label = 'Compton')
	plt.plot(ray_sum+comp_sum, label = 'Total scattering')
	plt.xlabel('collection semi-angle (deg)')
	plt.ylabel('cross section/(solid angle)$^2$ (cm$^2$/rad$^2$)')
	# plt.yscale('log')
	# plt.legend()
	# plt.ylim(1e-26,1e-22)
	plt.show()
	# plt.savefig('a.pdf',format='pdf')

if multiZ:
	ray_mat = []
	comp_mat = []
	for Z in range(1,101):			
		_dcs_rayleigh_pol = lambda beta_rad, theta_rad:dcs_rayleigh_pol(ev, Z, theta_rad, beta_rad)*np.sin(theta_rad)
		_dcs_compton_pol = lambda beta_rad, theta_rad:dcs_compton_pol(ev, Z, theta_rad, beta_rad)*np.sin(theta_rad)
		ray_sum = np.array([])
		comp_sum = np.array([])
		for x in range(1,91):
			omega=solid_angle(angle_range = [90-x, 90+x, -x, x], n_theta = 4*x, n_beta = 4*x)
			x = np.radians(x)
			a,b,gfun,hfun = [np.pi/2-x, np.pi/2+x, lambda y: -x, lambda y: x]
			ray = dblquad(_dcs_rayleigh_pol,a,b,gfun,hfun)[0]/omega.subtend**2
			comp = dblquad(_dcs_compton_pol,a,b,gfun,hfun)[0]/omega.subtend**2
			ray_sum = np.append(ray_sum,ray)
			comp_sum = np.append(comp_sum,comp)
		ray_mat.append(ray_sum)
		comp_mat.append(comp_sum)
	plt.figure(number2symbol(Z)+'_theta90')
	plt.title('cross section/(solid angle)$^2$ (cm$^2$/rad$^2$)')
	mat = np.array(ray_mat)+np.array(comp_mat)
	diff = mat[:,1:] - mat[:,0:-1]
	# plt.imshow(np.array(ray_mat)+np.array(comp_mat), norm = LogNorm())
	plt.imshow(diff)
	plt.xlabel('collection semi-angle (deg)')
	plt.ylabel('Z')
	plt.colorbar()
	plt.show()
		