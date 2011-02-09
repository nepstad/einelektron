import pyprop
import scipy
from pyprop.core import LmIndex
from numpy import conj, dot, abs, diff, r_, zeros, double, complex, array, linspace, pi
from numpy import maximum, sum, arctan2, imag, real, outer, sqrt, cos, sin, exp
from ..eigenvalues.eigenvalues import SetupRadialEigenstates, SetupOverlapMatrix
import eigenstates
from scipy.special import gamma, sph_harm
from scipy.interpolate import UnivariateSpline




def GetCoulombPhase(l, eta):
    """
    coulombPhase = GetCoulombPhase(l, eta)

    Generates the coloumb phase, arg[ gamma( l + 1 + i * eta) ].

    Parametres:
    -----------
    l : integer, orbital momentum quantum number.
    eta : double, Z/k. 
    """
    sigma = gamma(l + 1 + eta * 1j)
    return arctan2(sigma.imag, sigma.real)


class EigenstateAnalysis:
    """
    Perform analysis of a wavefunction in terms of eigenstates.
    """
    
    def __init__(self, conf):
	self.Config = conf
	self.BoundThreshold = 0.0
	self.m = 0
	    
    def Setup(self):
	#Setup pyprop problem.
	self.Problem = pyprop.Problem(self.Config)
	self.Problem.SetupStep()
	
	#Setup eigenstates/values.
	self.Eigenstate = eigenstates.Eigenstates(self.Config)
	
	#Setup overlap matrix.
	self.Overlap =  SetupOverlapMatrix(self.Problem)
	    
    def CalculateBoundDistribution(self, psi):
	"""
	Calculate norm**2 of projection onto all bound states. The given eigenvectors
	are assumed to cover all different angular momenta, but only for a given m (i.e. no 
	m-dependence in the potential).

	"""
	boundDistr = []
	boundE = []
	boundV = []
	boundTotal = 0
	
	for angIdx, curE, curV, l, m in self.Eigenstate.IterateBoundStates(self.BoundThreshold):
	    #Check for consistent psi
	    psiAngRange = psi.GetRepresentation().GetRepresentation(0).Range
	    assert psiAngRange.GetLmIndex(angIdx).l == l
	    assert psiAngRange.GetLmIndex(angIdx).m == m

	    #Get projection on eigenstates
	    psiSlice = psi.GetData()[angIdx, :]
	    overlapPsi = dot(self.Overlap, psiSlice)
	    proj = dot(conj(curV.transpose()), overlapPsi)

	    #Interpolate to get equispaced dP/dE
	    boundDistr.append( proj )
	    boundE.append( curE )
	    boundV.append(curV)
	    boundTotal += sum(abs(proj)**2)

	return boundE, boundV, boundDistr, boundTotal

    def CalculateBoundProbability(self, psi):
	"""
	boundTotal = CalculateBoundProbability(self, psi)

	Calculate norm squared of bound part of psi.

	Parametres
	----------
	psi : PyProp wavefunction object. The wavefunction after propagation.

	Returns
	-------
	boundTotal : double, the probability of being in a bound state.
	"""
	_, _, _, boundTotal = self.CalculateBoundDistribution(psi)
		
	return boundTotal
    
    

    def CalculateEnergyDistribution(self, psi):
	"""
	E, energyDistr = CalculateEnergyDistribution(self, psi)

	Calculates the energy distribution of the ejected electron.

	Parametres
	----------
	psi : PyProp wavefunction object. The wavefunction after propagation.

	Returns
	-------
	E : double 1D array, containing the energy grid.
	energyDistr : double 1D array, the energy distribution on E.
	"""
	#Defines energy grid.
	dE = min(diff(sorted(self.Eigenstate.EigenValues[0])))
	minE = 0
	maxE = self.Eigenstate.EigenValues[0][3*len(self.Eigenstate.EigenValues[0])/4]
	E = r_[minE:maxE:dE]
	
	#Initialise energy distribution array.
	energyDistr = zeros(len(E), dtype=double)
	
	#Initialise total ionisation.
	totalIon = 0

	#Loops over {l,m}, and corresponding eigenvalues/eigenvectors.
	
	for angIdx, curE, curV, l, m in self.Eigenstate.IterateStates(self.BoundThreshold):
	    #Get projection on eigenstates.
	    # | < V | S Psi > |**2
	    psiSlice = psi.GetData()[angIdx, :]
	    overlapPsi = dot(self.Overlap, psiSlice)
	    proj = abs(dot(conj(curV.transpose()), overlapPsi))**2
	    
	    #Add to total ionisation.
	    totalIon += sum(proj)

	    #Calculate density of states
	    interiorSpacing = list((diff(curE[:-1]) + diff(curE[1:])) / 2.)
	    #TODO:Lower dE at edges? Why?
	    leftSpacing = (curE[1] - curE[0]) / 2.
	    rightSpacing = (curE[-1] - curE[-2]) / 2.
	    spacing = array([leftSpacing] + interiorSpacing + [rightSpacing])
	    density = 1.0 / spacing

	    #Interpolate to get equispaced dP/dE
	    energyDistr +=  scipy.interp(E, curE, proj.ravel() * density, left=0, right=0) 
	
	#TODO:Controlling interpolation error. Remove when testd#|
	totalIon2 = sum(sum(array(energyDistr), axis=0)) * dE	#|
	print totalIon						#|
	print totalIon2						#|
	interpolateError = totalIon - totalIon2			#|

	return E, energyDistr


    def GetLegendrePoly(self, theta, phi):
	"""
	leg = GetLegendrePoly(lmIndices, theta)

	Calculates the Legendre polinomials for all {l,m}, and the angles in theta.

	Parametres
	----------
	lmIndices : list of lm-index objects in basis.
	theta : 1D double array, containing the theta grid.
	phi : 1D double array, containing the phi grid.

	Returns
	-------
	leg : 3D double array, containing legendre polinomials evaluated for different l, m and thetas.
	"""
	index_iterator = self.Config.AngularRepresentation.index_iterator

	leg = zeros([(index_iterator.lmax+1)**2, len(theta), len(phi)], dtype=complex)
	print "Legendre ..."
	for i, lm in enumerate(index_iterator.__iter__()):
	    print i
	    for j, my_theta in enumerate(theta):
		leg[i,j,:] = sph_harm(lm.m, lm.l, phi, my_theta)

	return leg


    def CalculateAngularDistribution(self, psi):
	"""
	theta, E, phi, angularDistribution = CalculateAngularDistribution(self, psi)

	Calculates the angular distribution.

	Parametres
	----------
	psi : PyProp wavefunction object. The wavefunction after propagation.

	Returns
	-------
	theta : 1D double array, containing the theta grid.
	E : 1D double array, containing the energy grid.
	phi : 1D double array, containing the phi grid.
	angularDistribution : 2D double array, containing the differentiated angular distribution.
	"""

	#Energy grid.
	#dE = maximum(min(diff(self.Eigenstate.EigenValues[0])), 0.05)
	dE = 0.01
	minE = dE
	#TODO: Insert intelligent value here.
	maxE = 5 #self.EigenValues[0][-1] #self.EigenValues[0][3*len(self.EigenValues[0])/4]
	E = r_[minE:maxE:dE]
	#Initialise energy distribution list.
	energyDistr = []
	
	#Charge.
	Z = 1
	
	#Angular grid.
	thetaCount = 100
	theta = linspace(0, pi, thetaCount)
	phiCount = 200
	phi = linspace(0, 2 * pi, phiCount)
	
	#Legendre polynomial values.
	leg = self.GetLegendrePoly(theta, phi)
	
	#Initialising the angular distribution array.
	angularDistrProj = zeros((thetaCount, len(E), phiCount), dtype=complex)
	    
	#Loops over ls, and corresponding eigenvalues/eigenvectors.
	#   curE is a 1D array,
	#   curV is a 2D array.	
	#Loops over {l,m}, and corresponding eigenvalues/eigenvectors.
	print "Projection ..."
	
	for angIdx, curE, curV, l, m in self.Eigenstate.IterateStates(self.BoundThreshold):
	    #From energy to momentum (k).
	    curk = sqrt(curE/2)

	    #Phase for outgoing waves
	    sigma = array([GetCoulombPhase(l, -Z / k) for k in curk])
	    phase = (-1.j)**l * exp(1.j * sigma)

	    #Get projection on eigenstates
	    psiSlice = psi.GetData()[angIdx, :]
	    overlapPsi = dot(self.Overlap, psiSlice)
	    proj = dot(conj(curV.transpose()), overlapPsi)

	    #Calculate density of states
	    #interiorSpacing = list((diff(curE[:-1]) + diff(curE[1:])) / 2.)
	    #TODO: Check versus EnergyDistribution.
	    interiorSpacing = list(diff(curE)[1:])
	    leftSpacing = (curE[1] - curE[0])
	    rightSpacing = (curE[-1] - curE[-2])
	    spacing = array([leftSpacing] + interiorSpacing + [rightSpacing])
	    density = 1.0 / sqrt(spacing)
	    
	    #The multiplying factors (except for the legendre pol.) of the distribution term. 
	    partialProj = phase * proj[:,0] * density 
		
	    #Interpolate (in complex polar coordinates) to get equispaced dP/dE.
	    r = abs(partialProj)**2
	    i = arctan2(imag(partialProj), real(partialProj))
	    argr = cos(i)
	    argi = sin(i)
	    
	    for ind in range(len(curE)):
		interpR    = UnivariateSpline(curE, r   , k=1, s=0)(E)
		interpArgR = UnivariateSpline(curE, argr, k=1, s=0)(E)
		interpArgI = UnivariateSpline(curE, argi, k=1, s=0)(E)

	    interpPhase = (interpArgR + 1.j*interpArgI) / sqrt(interpArgR**2 + interpArgI**2)
	    interpProj = sqrt(maximum(interpR, 0)) * interpPhase
	    
	    print angIdx


	    #Including the legendre pol., and adding the term to the complete distribution.
	    for ind in range(phiCount):
		angularDistrProj[:,:,ind] += outer(leg[angIdx,:,ind], interpProj)

	return theta, E, phi, abs(angularDistrProj)**2

