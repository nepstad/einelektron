"""
Returns the probability distribution of a one electron wavefunction in the xy-plane.
Main method is plot_wavefunction_xy.
"""
import sys
sys.path.append("../../")
sys.path.append("../../pyprop")
from scipy.interpolate import  RectBivariateSpline
from scipy.special import sph_harm

import einpartikkel
import pyprop
from numpy import linspace, pi, zeros, mod, outer, sqrt, arctan2

import einpartikkel.core.indexiterators
from einpartikkel.utils import UpdatePypropProjectNamespace
UpdatePypropProjectNamespace(pyprop.ProjectNamespace)


def plot_wavefunction_xy(filename,rmax):
	"""
	x_grid, y_grid, probability = plot_wavefunction_xy(filename, rmax)

	Returns the probability distribution of a one electron wavefunction in the xy-plane.

	Parametres
	----------
	filename : string, path and filename of the HDF% file containing the wavefunction.
	rmax : number, the radius to which you want to plot te wavefunction.

	Returns
	-------
	x_grid : 1D float array, the x values corresponding to 2nd axis of the probability distribution.
	y_grid : 1D float array, the y values corresponding to 1st axis of the probability distribution.
	probability : 2D float array, the probability distribution of a one electron wavefunction in the xy-plane.

	Example
	-------
	>>>import plot_wavefunction
	>>>from pylab import pcolormesh
	>>>x,y,P = plot_wavefunction.plot_wavefunction_xy("the/path/and/the/name_of_my_wavefunction.h5", 100)
	>>>pcolormesh(x,y,P)
	"""
	
	def grid_data(rmax,r_ratio,p_ratio):
		"""
		r_grid, phi_grid = grid_data(rmax,r_ratio,p_ratio)

		Define how the r & phi grid should be.

		Parametres
		----------
		rmax : number, size of the radial box for plotting.
		r_ratio : radial density of evaluation points.
		p_ratio : number of angular points will be rmax * p_ratio.

		Returns
		-------
		r_grid : the radial grid on which the wavefunction is evaluated.
		p_grid : the azimuthal grid (phi, in the xy-plane) on which the wavefunction is evaluated.
		"""
		r_grid = linspace(0,rmax,rmax*r_ratio)
		phi_grid = linspace(-1*pi,pi,rmax*p_ratio)
		return r_grid, phi_grid
	
	def convert_carthesian_to_polar(x,y):
		"""
		r,theta = convert_carthesian_to_polar(x,y)
		"""
		r = sqrt(x**2 + y**2)
		theta = arctan2(y,x)
		return r,theta

	#Load wavefunction and config object.
	psi = pyprop.CreateWavefunctionFromFile(filename)#,datasetPath="/initialWavefunction")
	psi_data = psi.GetData()
	conf = pyprop.LoadConfigFromFile(filename)
	
	#Create B-spline object.
	spline_object = pyprop.BSPLINE()
	spline_object.ApplyConfigSection(conf.RadialRepresentation)
	
	#Get grid data.
	r_ratio = 2
	p_ratio = 3
	r_grid, phi_grid = grid_data(rmax,r_ratio,p_ratio)
	
	
	#Evaluate spherical harmonics.
	leg = GetLegendrePoly(conf, [pi/2.], phi_grid)
	
	#Initialise grid for wavefunction.
	psi_grid = zeros([len(r_grid), len(phi_grid)], dtype = "complex")
	
	#Wavefunction to grid.
	for i in range(leg.shape[0]):
		temp_rad = spline_object.ConstructFunctionFromBSplineExpansion(psi_data[i,:], r_grid)
		psi_grid += outer(temp_rad,leg[i, 0,:])
	
	prob_grid = abs(psi_grid)**2

	#Interpolation function.
	interp_prob = RectBivariateSpline(r_grid, phi_grid, prob_grid, kx=1, ky=1, s = 0)
	
	xmax = rmax/sqrt(2)
	#To get a denser grid in xy, increase xsize.
	xsize = xmax * 2
	x_grid = linspace(-1 * xmax, xmax, xsize)
	
	#To avoid mixup this is different from x. After debug, feel free to change this.
	y_grid = linspace(-1 * xmax, xmax, xsize)
	
	#Notice: y on first axis (because of plotting).
	prob_xgrid = zeros([len(y_grid), len(x_grid)])

	for ind_x, x in enumerate(x_grid):
		for ind_y, y in enumerate(y_grid):
			r, theta = convert_carthesian_to_polar(x, y)
			prob_xgrid[ind_y, ind_x] = interp_prob(r, theta)
	
	#Interpolate to denser grid.	
	x_grid_dense = linspace(-1 * xmax, xmax, 5*xsize)
	y_grid_dense = linspace(-1 * xmax, xmax, 5*xsize)
	prob_dense = RectBivariateSpline(x_grid, y_grid, prob_xgrid, kx=1, ky=1, s = 0)(x_grid_dense, y_grid_dense)

	return x_grid_dense, y_grid_dense, prob_dense
			


def GetLegendrePoly(conf, theta, phi):
	"""
	leg = GetLegendrePoly(conf, theta, phi)

	Calculates the Legendre polynomials for all {l,m}, and the angles in theta.

	Parametres
	----------
	conf : Config object.
	theta : 1D double array, containing the theta grid.
	phi : 1D double array, containing the phi grid.

	Returns
	-------
	leg : 3D double array, containing legendre polinomials evaluated for different l, m and thetas.
	"""
	index_iterator = conf.AngularRepresentation.index_iterator
	nr_lm = (index_iterator.lmax+1)**2

	leg = zeros([nr_lm, len(theta), len(phi)], dtype=complex)
	print "Legendre ..."
	for i, lm in enumerate(index_iterator.__iter__()):
		if mod(i,nr_lm/10)==0:
			print i * 100. / nr_lm, "%"
		for j, my_theta in enumerate(theta):
			leg[i,j,:] = sph_harm(lm.m, lm.l, phi, my_theta)

	return leg
