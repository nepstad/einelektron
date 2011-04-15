import pyprop
import scipy
from numpy import conj, dot, abs, diff, r_, zeros, double, complex, array, linspace, pi
from numpy import maximum, sum, arctan2, imag, real, outer, sqrt, cos, sin, exp, shape
import eigenstates
import coulombwaves
import os.path
import tables
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
		
		#Multiply overlap on wavefunction
		overlapPsi = self.MultiplyOverlap(psi)
		
		for angIdx, curE, curV, l, m in self.Eigenstate.IterateBoundStates(self.BoundThreshold):
			#Get projection on eigenstates
			psiSlice = overlapPsi.GetData()[angIdx, :]
			proj = dot(conj(curV.transpose()), psiSlice)

			#Interpolate to get equispaced dP/dE
			boundDistr.append( proj )
			boundE.append( curE )
			boundV.append(curV)
			boundTotal += sum(abs(proj)**2)

		return boundE, boundV, boundDistr, boundTotal


	def CalculateAllBoundProjections(self, psi):
		"""
		Calculate norm**2 of projection onto all bound states. Returns the 
		projection on each state, specified by energy and {l,m}.
		(In the case of hydrogen the energy is easily mappet to the n-quantum number.)


		Parametres
		----------
		psi : PyProp wavefunction object. The wavefunction after propagation.

		Returns
		-------
		boundDistr : list, the projections onto all boundstate.
		lmIndices : list, the {l,m}'s of of each boundstate.
		energy : list, the energy of each boundstate.

		"""
		boundDistr = []
		lmIndices = []
		energy = []
		
		#Multiply overlap on wavefunction
		overlapPsi = self.MultiplyOverlap(psi)
		
		for angIdx, curE, curV, l, m in self.Eigenstate.IterateBoundStates(self.BoundThreshold):
			#Get projection on eigenstates
			psiSlice = overlapPsi.GetData()[angIdx, :]
			proj = dot(conj(curV.transpose()), psiSlice)


			#Stuff eigenvalues, projections and {l,m} into a list
			for E,p in zip(curE, proj):
				boundDistr.append(abs(p)**2)
				lmIndices.append((l,m))
				energy.append(E)


		return boundDistr, lmIndices, energy


	def RemoveProjectionOnBoundStates(self, psi):
		"""Remove projection on bound states in-place

		"""

		#Multiply overlap on wavefunction
		overlapPsi = self.MultiplyOverlap(psi)

		for angIdx, curE, curV, l, m in self.Eigenstate.IterateBoundStates(self.BoundThreshold):
			#Get projection on eigenstates
			psiSlice = overlapPsi.GetData()[angIdx, :]
			proj = dot(conj(curV.transpose()), psiSlice)
			for p, v in zip(proj, curV.transpose()):
				psi.GetData()[angIdx, :] -= p * v[0,:]


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
    

	def CalculateEnergyDistribution(self, psi, minE, maxE, dE):
		"""
		E, energyDistr = CalculateEnergyDistribution(self, psi)

		Calculates the energy distribution of the ejected electron.

		Parametres
		----------
		psi:  PyProp wavefunction object. The wavefunction after propagation.
		minE: lower cutoff in energy spectrum
		maxE: upper cutoff in energy spectrum
		dE:   spacing of interpolated energies in spectrum

		Returns
		-------
		E : double 1D array, containing the energy grid.
		energyDistr : double 1D array, the energy distribution on E.

		"""

		self.Logger.info("Now calculating dP/dE...")

		#Multiply overlap on wavefunction
		overlapPsi = self.MultiplyOverlap(psi)

		#Create energy grid.
		E = r_[minE:maxE:dE]
		
		#Initialise energy distribution array.
		energyDistr = zeros(len(E), dtype=double)
		
		#Initialise total ionisation.
		totalIon = 0

		#Loops over {l,m}, and corresponding eigenvalues/eigenvectors.
		for angIdx, curE, curV, l, m in self.Eigenstate.IterateStates(self.BoundThreshold):
			#Get projection on eigenstates.
			# | < V | S Psi > |**2
			psiSlice = overlapPsi.GetData()[angIdx, :]
			proj = abs(dot(conj(curV.transpose()), psiSlice))**2
			
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



	def CreateAndSave_l_matrices(self, lmin, filename, theta, phi):
		"""
		CreateAndSave_l_matrices(lmin, filename, theta, phi)
		
		Calculates the Legendre polynomials for all {l,m} and all angles.
		The polynomials are stored in a h5-file. If the file exists, the
		file is updated.

		Parametres
		----------
		lmin : int, calculate Legendre polynomials for l >= lmin.
		filename : string, grid file name.
		theta : 1D double array, containing the theta grid.
		phi : 1D double array, containing the phi grid.

		"""
		if lmin == 0:
			mode = 'w'
		else:
			mode = 'r+'

		f = tables.openFile(filename, mode)
		root = f.root
		index_iterator = self.Config.AngularRepresentation.index_iterator

		print "Legendre ..."
		prevl = -1;
		for i, lm in enumerate(index_iterator.__iter__()):
			print i
			if lm.l >= lmin:
				if lm.l != prevl:
					midx = 0
					leg = zeros([(2 * lm.l + 1), len(theta), len(phi)], dtype=complex)

				for j, my_theta in enumerate(theta):
					leg[midx,j,:] = sph_harm(lm.m, lm.l, phi, my_theta)
				
				midx += 1

				if midx == 2 * lm.l + 1:
					f.createArray('/','l_' + str(lm.l),leg)

				prevl = lm.l
		f.setNodeAttr("/","lmax",index_iterator.lmax)
		f.close()

	
	def LoadAndCreateLegendreMatrix(self, filename, theta, phi):
		"""
		leg = LoadAndCreateLegendreMatrix(filename, theta, phi)
		
		Loads data from a file created by CreateAndSave_l_matrices(), and creates
        a 3D double array, containing legendre polinomials evaluated for different l, m and thetas.

		Parametres
		----------
		filename : string, grid file name.
		theta : 1D double array, containing the theta grid.
		phi : 1D double array, containing the phi grid.

		Returns
		-------
		leg : 3D double array, containing legendre polinomials evaluated for different l, m and thetas.

		"""

		print "Load Legendre polys..."
		index_iterator = self.Config.AngularRepresentation.index_iterator

		leg = zeros([(index_iterator.lmax+1)**2, len(theta), len(phi)], dtype=complex)	

		lmIdx = 0
		f = tables.openFile(filename)
		for l in range(0, index_iterator.lmax + 1):
			a = eval('f.root.' + 'l_' + str(l))
			shp = shape(a)
			leg[lmIdx:lmIdx + shp[0],:,:] = a[:]

			lmIdx += shp[0]

		f.close()

		return leg


	def GetLegendrePoly(self, theta, phi):
		"""
		leg = GetLegendrePoly(lmIndices, theta)

		Calculates the Legendre polinomials for all {l,m}, and the angles in theta.
		The polynomials are stored in a grid spesific file. If the file already exists, the file is updated.

		Parametres
		----------
		lmIndices : list of lm-index objects in basis.
		theta : 1D double array, containing the theta grid.
		phi : 1D double array, containing the phi grid.

		Returns
		-------
		leg : 3D double array, containing legendre polinomials evaluated for different l, m and thetas.
		"""

		gridDir = 'legendre'
		gridFile = gridDir + '/SphericalGrid_%.4f_%.4f_%i_%.4f_%.4f_%i.h5'%(theta[0],theta[-1],len(theta),phi[0],phi[-1],len(phi))

		index_iterator = self.Config.AngularRepresentation.index_iterator

		#Check if the grid already exists, if not, create a file.
		if os.path.exists(gridFile):

			f = tables.openFile(gridFile)
			lmaxStored = f.root._v_attrs.lmax
			f.close()

			#Update?
			if index_iterator.lmax > lmaxStored:
				print "The file " + gridFile + " is updated to support lmax " + str(index_iterator.lmax)
				self.CreateAndSave_l_matrices(lmaxStored + 1, gridFile, theta, phi)	

		else:
			if not os.path.exists(gridDir):
				os.mkdir(gridDir)

			self.CreateAndSave_l_matrices(0, gridFile, theta, phi)


		#Load the legendre polys from file
		leg = self.LoadAndCreateLegendreMatrix(gridFile, theta, phi)

		return leg


#	def GetLegendrePoly(self, theta, phi):
#		"""
#		leg = GetLegendrePoly(lmIndices, theta)
#
#		Calculates the Legendre polinomials for all {l,m}, and the angles in theta.
#
#		Parametres
#		----------
#		lmIndices : list of lm-index objects in basis.
#		theta : 1D double array, containing the theta grid.
#		phi : 1D double array, containing the phi grid.
#
#		Returns
#		-------
#		leg : 3D double array, containing legendre polinomials evaluated for different l, m and thetas.
#		"""
#		index_iterator = self.Config.AngularRepresentation.index_iterator

#		leg = zeros([(index_iterator.lmax+1)**2, len(theta), len(phi)], dtype=complex)
#		print "Legendre ..."
#		for i, lm in enumerate(index_iterator.__iter__()):
#			print i
#			for j, my_theta in enumerate(theta):
#				leg[i,j,:] = sph_harm(lm.m, lm.l, phi, my_theta)
#
#		return leg




	def CalculateAngularDistribution(self, psi, minE, maxE, dE, Z=1):
		"""
		theta, E, phi, angularDistribution = CalculateAngularDistribution(self, psi)

		Calculates the angular distribution.

		Parametres
		----------
		psi : PyProp wavefunction object. The wavefunction after propagation.
		Z:    (int) Coulomb wave charge

		Returns
		-------
		theta : 1D double array, containing the theta grid.
		E : 1D double array, containing the energy grid.
		phi : 1D double array, containing the phi grid.
		angularDistribution : 2D double array, containing the differentiated angular distribution.

		"""

		#Multiply overlap on wavefunction
		overlapPsi = self.MultiplyOverlap(psi)

		E = r_[minE:maxE:dE]
		#Initialise energy distribution list.
		energyDistr = []
		
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
			curk = sqrt(curE*2.)

			#Phase for outgoing waves
			sigma = array([GetCoulombPhase(l, -Z / k) for k in curk])
			phase = (-1.j)**l * exp(1.j * sigma)

			#Get projection on eigenstates
			psiSlice = overlapPsi.GetData()[angIdx, :]
			
			proj = dot(conj(curV.transpose()), psiSlice)

			#Calculate density of states
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

#			#TODO REMOVE
#			return proj, angIdx, curE, outer(leg[angIdx,:,ind], interpProj)
#			#

		return theta, E, phi, abs(angularDistrProj)**2


	def MultiplyOverlap(self, inPsi):
		"""Multiply overlap matrix on wavefunction, return new wavefunction

		"""
		#Get a copy of psi
		outPsi = inPsi.Copy()

		#Muliply overlap
		outPsi.GetRepresentation().MultiplyOverlap(outPsi)

		return outPsi


class CoulombwaveAnalysis(EigenstateAnalysis):
	"""
	Perform analysis of a wavefunction in terms of Coulomb waves.
	"""

	def __init__(self, conf, Z, fileList):
		self.Config = conf
		self.BoundThreshold = 0.0
		self.m = 0
		self.Z = Z
		self.FileList = fileList
	    
	def Setup(self):
		#Setup Coulomb waves
		#folderName = "CoulombWaves"
		#fname = coulombwaves.CoulombWaveNameGenerator(self.Config, folderName,\
		#		self.Z, dE, Emax)
		self.Eigenstate = coulombwaves.CoulombWaves.FromMultipleFiles(self.FileList)
		
		#Check that the Coulomb waves are compatible (TODO: more checks)
		assert self.Eigenstate.Z == self.Z

