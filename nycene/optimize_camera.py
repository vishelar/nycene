from scipy.optimize import minimize
import numpy as np
from random import gauss
from colin import colin
import fiducials



def fullfunc(params,xyz_s,xy_t):

	''' Find the sum of squares difference '''
	omega, phi, kappa, xs, ys, zs, f = params

#	if (omega<0.0) or (omega>=2.0*np.pi):
	if (omega<0.0) or (omega>=2.0*np.pi):
		return 1e9+omega**4
	elif (phi<-0.5*np.pi) or (phi>=0.5*np.pi):
		return 1e9+phi**4
	elif (kappa<0.0) or (kappa>=2.0*np.pi):
		return 1e9+kappa**4
	elif zs<0.0:
		return 1e9+zs**4
	elif f<0.0:
		return 1e9+f**4
	elif (np.abs(params[3] - 977119)>1000) or \
		    (np.abs(params[4] - 210445)>1000):
		return 1e9 + xs**2

#	params[3] = 977119.
#	params[4] = 210445.

#	colin_xy = -1.0*colin(params,xyz_s.astype(float))
	colin_xy = 1.0*colin(params,xyz_s.astype(float))
	diff = ((colin_xy - xy_t)**2).sum()

	return diff



def random_start(params):

	''' Perturbs camera position in a gaussian fashion '''
	return params
#	return params + np.array([gauss(0, 0.1), gauss(0, 0.1), \
##	gauss(0, 0.1),gauss(0, 20000),gauss(0, 20000), \
#	gauss(0, 0.1),gauss(0, 100),gauss(0, 100), \
#	gauss(0, 100), gauss(0, 30)])
##	gauss(0, 10), gauss(0, 2000)])



def call(params,xyz_s,xy_t):

	''' Guess parameters near start and brute-force minimize '''
	start = random_start(params)
#	print ("start, xyz_s, xy_t",start,xyz_s,xy_t)
	res = minimize(fullfunc, start,args=(xyz_s,xy_t), \
#		method = 'Powell', \
#		options={'ftol':1e-6})
		method = 'Nelder-Mead', \
		options={'maxfev': 10000, 'maxiter': 10000})

	return res



def run(imname, num_iter, params=None):

	''' Run the optimization routine for the given image '''

	if params == None:
		print "Choosing 1MT UO as initial guess\n"
		params = np.array([1.56926535e+00, -1.20789690e-01,\
			-3.05255789e-03, 9.87920425e+05, 1.91912958e+05, \
			3.85333237e+02, -1.10001068e+04]) 

	xyz_s = fiducials.lidar_fiducials(imname)
	xy_t = fiducials.image_fiducials(imname,center=True)

	min_score = 100000000000000

	for i in range(0, num_iter):
		result = call(params,xyz_s,xy_t)
		print "params, score", result.x, result.fun
#		print('doing it {0} {1}'.format(i,result.fun))
#		import pdb; pdb.set_trace()
		if (result.fun < min_score):# and (result.x[3] < 980491):
			min_score = result.fun
			params = result.x

	return [min_score, params]



