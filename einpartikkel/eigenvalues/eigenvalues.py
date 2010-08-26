from numpy import *
import scipy
import scipy.linalg

import pyprop

from ..utils import RegisterAll

@RegisterAll
def SetupRadialEigenstates(prop, potentialIndices=[0], mList = [0]):
	"""
	Finds the eigenvalues and eigenvectors of the given potentials
	of prop. 

	The eigenvalues are found by setting up a radial matrix for each (l,m)-value
	and using the generalized eigenvalue solver in scipy to find all
	eigenvalues and vectors. 
	
	eigenvalues is a list of 1-d eigenvalue arrays. 
	
	"""

	if not pyprop.IsSingleProc():
		raise Exception("Works only on a single processor")

	S = SetupOverlapMatrix(prop)
	bspl = prop.psi.GetRepresentation().GetRepresentation(1).GetBSplineObject()
	angRepr = prop.psi.GetRepresentation().GetRepresentation(0)
	phaseGrid = array((0, bspl.GetBreakpointSequence()[1]), dtype=double)
	phaseBuffer = zeros(2, dtype=complex)

	eigenValues = []
	eigenVectors = []
	angIdxList = []
	lmIdxList = []

	lmCount = prop.psi.GetData().shape[0]
	eigenvectorScaling = 1

	for angIdx in range(lmCount):
		lmIdx = angRepr.Range.GetLmIndex(angIdx)
		if not lmIdx.m in mList:
			continue
		print "Calculating eigenpairs for l,m = %s,%s..." % (lmIdx.l, lmIdx.m)
		M = SetupRadialMatrix(prop, potentialIndices, angIdx)

		E, V = scipy.linalg.eig(a=M, b=S)

		idx = argsort(real(E))
		E = real(E[idx])
		eigenValues.append(E)

		#Sort and normalize eigenvectors
		sNorm = lambda v: sqrt(abs(sum(conj(v) * dot(S, v))))
		V = array([v/sNorm(v) for v in [V[:,idx[i]] for i in range(V.shape[1])]]).transpose()
		eigenVectors.append(V)

		#assure correct phase convention (first oscillation should start out real positive)
		for i, curE in enumerate(E):
			bspl.ConstructFunctionFromBSplineExpansion(V[:,i].copy(), phaseGrid, phaseBuffer)
			phase = arctan2(imag(phaseBuffer[1]), real(phaseBuffer[1]))
			V[:,i] *= exp(-1.0j * phase)

		angIdxList += [angIdx]
		lmIdxList += [lmIdx]

	return eigenValues, eigenVectors, angIdxList, lmIdxList

@RegisterAll
def SetRadialEigenstate(psi, eigenVectors, angIdxList, quantumNumbers, sourceScaling=0., destScaling=1.0):
	"""
	Sets psi to an eigenvector from a list of eigenvectors as calculated by
	SetupRadialEigenstates()

	"""
	lmIdx = quantumNumbers.GetLmIndex()
	angIdx = psi.GetRepresentation().GetRepresentation(0).Range.GetGridIndex(lmIdx)
	if not angIdx in angIdxList:
		raise Exception("That eigenstate is not included in 'eigenVectors!")
	radialIndex = quantumNumbers.GetRadialIndex()
	eigAngIdx = angIdxList.index(angIdx)
	vec = eigenVectors[eigAngIdx][:, radialIndex]
	psi.GetData()[:] *= sourceScaling
	psi.GetData()[angIdx, :] += destScaling * vec


def SetupRadialMatrix(prop, whichPotentials, angularIndex):
	matrixSize = prop.psi.GetData().shape[1]
	matrix = zeros((matrixSize, matrixSize), dtype=complex)

	for potNum in whichPotentials:	
		if isinstance(potNum, pyprop.TensorPotential):
			potential = potNum
		else:
			potential = prop.Propagator.BasePropagator.PotentialList[potNum]
		print "    Processing potential: %s" % (potential.Name, )

		angularBasisPairs = potential.BasisPairs[0]
		idx = [idx for idx, (i,j) in enumerate(zip(angularBasisPairs[:,0], angularBasisPairs[:,1])) if i==j==angularIndex]
		if len(idx) != 1:
			raise "Invalid angular indices %s" % idx
		idx = idx[0]

		basisPairs = potential.BasisPairs[1]

		for i, (x,xp) in enumerate(basisPairs):
			indexLeft = x
			indexRight = xp
			matrix[indexLeft, indexRight] += potential.PotentialData[idx, i]

	return matrix


def SetupOverlapMatrix(prop):
	overlap = prop.Propagator.BasePropagator.GeneratePotential(prop.Config.OverlapPotential)
	overlap.SetupStep(0.)
	matrix = SetupRadialMatrix(prop, [overlap], 0)
	return matrix
