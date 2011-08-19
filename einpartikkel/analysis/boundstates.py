"""
boundstates
===========

Basic bound state operations are provided. These are:

	1)Load bound states from disk.
	2)Remove projection on bound states

"""
import numpy
from numpy import zeros, complex, array, dot, conj
import pyprop
from ..utils import RegisterAll
from ..eigenvalues.eigenvalues_iter import LoadEigenpairs
from ..eigenvalues.eigenvalues import SetupRadialEigenstates, SetupOverlapMatrix


@RegisterAll
class Boundstates(object):
	"""A simple class to represent bound states and perform project-remove

	"""

	def __init__(self, conf):
		self.Config = conf
		self.LmList = [(lmIdx.l, lmIdx.m) for lmIdx in \
				conf.AngularRepresentation.index_iterator]
		self.Energies = []
		self.States = []
		self.Psi = None
		self.Overlap = None
		self.IsSetup = False

		self.Logger = pyprop.GetClassLogger(self)

	def Setup(self):
		"""Load eigenstates from disk
		"""

		#Load states
		for l,m in self.LmList:
			curE, curV = LoadEigenpairs(self.Config, l, m)
			if len(curE) > 0:
				idx = numpy.where(curE < 0.0)[0]
				curE = array(curE)[idx]
				curV = array(curV)[idx]
			self.Energies += [curE]
			self.States += [curV]

		self.Psi = pyprop.CreateWavefunction(self.Config)

		#Setup overlap matrix
		#self.Overlap = SetupOverlapMatrix(self.Config.OverlapPotential, \
		#		self.Psi)

		self.IsSetup = True

	def RemoveProjection(self, psi):
		"""Remove projection of psi onto all bound states

		"""
		assert self.IsSetup

		#Copy psi
		psiTmp = psi.Copy()

		#Multiply overlap
		#radialRank = 1
		psiTmp.GetRepresentation().MultiplyOverlap(psiTmp)

		angRange = psi.GetRepresentation().GetRepresentation(0).Range
		for idx, (l,m) in enumerate(self.LmList):
			if len(self.States) == 0:
				continue
			self.Logger.info("Removing projection on (l,m) = (%i,%i)" % (l,m))
			angIdx = angRange.GetGridIndex(pyprop.core.LmIndex(l,m))
			assert (angIdx > -1)
			curV = self.States[idx]
			if len(curV) == 0:
				continue
			#psiSlice = psi.GetData()[angIdx, :]
			#overlapPsi = dot(self.Overlap, psiSlice)
			#proj = dot(conj(curV), overlapPsi)
			psiSlice = psiTmp.GetData()[angIdx, :]
			proj = dot(conj(curV), psiSlice)
			#print "l=%s, m=%s, p=%s" % (l,m, proj)
			for p, v in zip(proj, curV):
				psi.GetData()[angIdx, :] -= p * v[:]


def SetupOverlapMatrix(confSec, psi):
	"""Setup overlap matrix (full)

	"""
	#Create TensorPotential
	overlap = pyprop.TensorPotential(psi)
	confSec.Apply(overlap)
	overlap.SetupStep(0.)

	#Transform to full matrix representation
	matrix = SetupRadialMatrix(psi, overlap, 0)

	return matrix


def SetupRadialMatrix(psi, potential, angularIndex):
	"""Create standard (full) matrix from a tensor potential
	"""
	matrixSize = psi.GetData().shape[1]
	matrix = zeros((matrixSize, matrixSize), dtype=complex)

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
