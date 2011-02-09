"""
eigenvalues
===========

Eigenvalue calculations module. Provides high-level interface to iterative
methods for finding some eigenvalues of a large Hamiltonian


Special config section
----------------------
A special config section is required by this module, to facilitate computations
and serialization of eigenpairs:

[Eigenvalues]
folder_name = "/path/to/eigenpairs  #where to store eigenpairs, unique to
									#potential 
shift = -0.5                        #shift used in shift-invert procedure

"""

import os
import sys
from numpy import where, array
import tables
import pypar
import pyprop
from pyprop import GetFunctionLogger
from pyprop.serialization import RemoveExistingDataset
from einpartikkel.utils import RegisterAll 
from einpartikkel.namegenerator import GetRadialPostfix

@RegisterAll
def FindEigenvaluesInverseIterationsPiram(conf):
	"""
	Calculates eigenvalues and eigenvectors for conf around a shift
	by using inverse iterations and pIRAM.

	Input
	-----
	conf: a pyprop config object

	Returns
	-------
	Piram solver object
	GMRES solver object
	Computed eigenvalues

	"""

	#Setup Problem
	prop = pyprop.Problem(conf)
	prop.SetupStep()

	#Setup shift invert solver in order to perform inverse iterations
	shiftInvertSolver = pyprop.GMRESShiftInvertSolver(prop)
	prop.Config.Arpack.inverse_iterations = True
	prop.Config.Arpack.matrix_vector_func = shiftInvertSolver.InverseIterations

	#Setup eiganvalue solver
	solver = pyprop.PiramSolver(prop)
	solver.Solve()

	#Get the converged eigenvalues
	ev = solver.Solver.GetEigenvalues().copy()
	estimates = solver.Solver.GetConvergenceEstimates().copy()
	idx = where(estimates < 0)[0]
	eigenvalues = ev[idx]

	#convert from shift inverted eigenvalues to "actual" eigenvalues
	shift = conf.Eigenvalues.shift
	eigenvalues = 1.0 / eigenvalues + shift

	return solver, shiftInvertSolver, eigenvalues


@RegisterAll
def SaveEigenvalueSolverShiftInvert(solver, shiftInvertSolver):
	"""
	Saves the output of FindEigenvaluesNearShift, including error estimates 
	to a hdf file.
	"""
	
	logger = GetFunctionLogger()

	conf = solver.BaseProblem.Config
	shift = conf.Eigenvalues.shift
	angularRank = 0

	#Get angular momentum and z-projection
	#l = conf.AngularRepresentation.index_iterator.l
	#m = conf.AngularRepresentation.index_iterator.m
	psi = solver.BaseProblem.psi
	angRange = \
		psi.GetRepresentation().GetRepresentation(angularRank).Range
	assert angRange.GetGrid().size == 1
	l = angRange.GetLmIndex(0).l
	m = angRange.GetLmIndex(0).m

	#generate filename
	filename = NameGeneratorBoundstates(conf, l, m)

	#Get eigenvalue error estimates
	errorEstimatesPIRAM = solver.Solver.GetErrorEstimates()
	convergenceEstimatesEig = solver.Solver.GetConvergenceEstimates()
	errorEstimatesGMRES = shiftInvertSolver.Solver.GetErrorEstimateList()

	#Get eigenvalues
	prop = solver.BaseProblem
	E = 1.0 / array(solver.GetEigenvalues()) + shift

	#remove file if it exists
	try:
		if os.path.exists(filename):
			if pyprop.ProcId == 0:
				os.remove(filename)
	except:
		logger.error("Could not remove %s (%s)" % (filename, sys.exc_info()[1]))

	#Store eigenvalues and eigenvectors
	logger.info("Now storing eigenvectors...")
	for i in range(len(E)):
		solver.SetEigenvector(prop.psi, i)
		prop.SaveWavefunctionHDF(filename, GetEigenvectorDatasetPath(i))

	if pyprop.ProcId == 0:
		eigGroupName = GetEigenvectorGroupName()
		RemoveExistingDataset(filename, "/%s/Eigenvalues" % eigGroupName)
		RemoveExistingDataset(filename, "/%s/ErrorEstimateListGMRES" \
				% eigGroupName)
		RemoveExistingDataset(filename, "/%s/ErrorEstimateListPIRAM" \
				% eigGroupName)
		RemoveExistingDataset(filename, "/%s/ConvergenceEstimateEig" \
				% eigGroupName)
		h5file = tables.openFile(filename, "r+")
		try:
			#myGroup = h5file.createGroup("/", "Eig")
			myGroup = h5file.getNode(eigGroupName)
			h5file.createArray(myGroup, "Eigenvalues", E)
			h5file.createArray(myGroup, "ErrorEstimateListGMRES", errorEstimatesGMRES)
			h5file.createArray(myGroup, "ErrorEstimateListPIRAM", errorEstimatesPIRAM)
			h5file.createArray(myGroup, "ConvergenceEstimateEig", convergenceEstimatesEig)

			#Store config
			myGroup._v_attrs.configObject = prop.Config.cfgObj
			
			#PIRAM stats
			myGroup._v_attrs.opCount = solver.Solver.GetOperatorCount()
			myGroup._v_attrs.restartCount = solver.Solver.GetRestartCount()
			myGroup._v_attrs.orthCount = solver.Solver.GetOrthogonalizationCount()
		except:
			logger.warning("Warning: could not store eigenvalues and error estimates!")
		finally:
			h5file.close()

		pypar.barrier()


@RegisterAll
def NameGeneratorBoundstates(conf, l, m):
	"""Returns a file name generated from a Pyprop config 

	Parametres
	----------
	conf: config object.
	l:    (int) angular momentum
	m:    (int) m quantum number

	Returns
	-------
	fileName : string, excellent file name.

	"""
	#Folder
	folderName = conf.Eigenvalues.folder_name
	
	#Grid characteristics
	radialPostfix = "_".join(GetRadialPostfix(conf))
	
	filename = "%s/" % folderName
	filename += radialPostfix + "_l%i_m%i" % (l,m) + ".h5"
	
	return os.path.normpath(filename)


@RegisterAll
def GetEigenvectorDatasetName(idx):
	return "Eigenvector_%i" % idx

@RegisterAll
def GetEigenvectorGroupName():
	return "/Eig"

@RegisterAll
def GetEigenvectorDatasetPath(eigenvectorIndex):
	return "%s/%s" % (GetEigenvectorGroupName(), \
		GetEigenvectorDatasetName(eigenvectorIndex))

