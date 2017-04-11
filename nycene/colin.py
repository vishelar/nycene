import numpy as np
import scipy as sp
from scipy.weave import inline



def colin(params, xyz_a):

	# Unwrap params
	kappa, phi, omega, xs, ys, zs, f = params

	omega = float(omega)
	phi = float(phi) + 0.5*np.pi
	kappa = float(kappa)
	xs = float(xs)
	ys = float(ys)
	zs = float(zs)
	f = float(f)

	# -- utils
	co = np.cos(omega)
	so = np.sin(omega)
	cp = np.cos(phi)
	sp = np.sin(phi)
	ck = np.cos(kappa)
	sk = np.sin(kappa)

	a1 =  cp*ck+sp*so*sk
	b1 =  cp*sk+sp*so*ck
	c1 =  sp*co
	a2 = -co*sk
	b2 =  co*ck
	c2 =  so
	a3 =  sp*ck+cp*so*sk
	b3 =  sp*sk-cp*so*ck
	c3 =  cp*co

	ynum  = a1*(xyz_a[:,0]-xs)+b1*(xyz_a[:,1]-ys)+c1*(xyz_a[:,2]-zs)
	xnum  = a2*(xyz_a[:,0]-xs)+b2*(xyz_a[:,1]-ys)+c2*(xyz_a[:,2]-zs)
	denom = a3*(xyz_a[:,0]-xs)+b3*(xyz_a[:,1]-ys)+c3*(xyz_a[:,2]-zs)

	xx = -f*xnum/denom
	yy = f*ynum/denom

	return np.vstack([xx,yy]).T

        '''	
	# Get number of fiducial points
	N = int(xyz_a.shape[0])

	# Initialize the result vector
	colin_xy = np.zeros((N,2))

	# Define c code that will evaluate the functions
	code = """

		double sinp = sinf(phi);
		double cosp = cosf(phi);
		double sino = sinf(omega);
		double coso = cosf(omega);
		double sink = sinf(kappa);
		double cosk = cosf(kappa);

		double a1 = cosp*cosk;
		double b1 = coso*sink + sino*sinp*cosk;
		double c1 = sino*sink - coso*sinp*cosk;
		double a2 = -1 * cosp*sink;
		double b2 = coso*cosk - sino*sinp*sink;
		double c2 = sino*cosk + coso*sinp*sink;
		double a3 = sinp;
		double b3 = -1*sino*cosp;
		double c3 = cosp*coso;

		int j;
		int k;
		double denom;

		for(int i = 0; i < N; i++)
		{
			j = i*2;
			k = i*3;
			denom = (a3*(xyz_a[k]-xs) + b3*(xyz_a[(k+1)]-ys) + c3*(xyz_a[(k+2)]-zs));

			colin_xy[j] = -1.0*f*(a1*(xyz_a[k]-xs) + b1*(xyz_a[(k+1)]-ys) + c1*(xyz_a[(k+2)]-zs))/denom;
			colin_xy[(j+1)] = -1.0*f*(a2*(xyz_a[k]-xs) + b2*(xyz_a[(k+1)]-ys) + c2*(xyz_a[(k+2)]-zs))/denom;
		}
		return_val = 1;
	"""

	# Use scipy.weave.inline to run the c code
	res = inline(code, ['colin_xy', 'omega', 'phi', 'kappa', 'xs', \
		'ys', 'zs', 'f', 'xyz_a', 'N'], headers = ['<math.h>'], \
		compiler = 'gcc')
	
	# Return the pixel (x,y) positions
	return colin_xy
	'''
