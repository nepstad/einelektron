import pyprop
from pyprop.core import LmIndex
from numpy import dot, conj, where
from ..eigenvalues.eigenvalues import SetupRadialEigenstates, SetupOverlapMatrix


def CalculateBoundDistribution(psi, eigenValues, eigenVectors, lmIdxList, m, overlap, boundThreshold=0):
	"""
	Calculate norm**2 of projection onto all bound states. The given eigenvectors
	are assumed to cover all different angular momenta, but only for a given m (i.e. no 
	m-dependence in the potential).

	"""

	angularRank = 0

	boundDistr = []
	boundE = []
	boundV = []
	boundTotal = 0

	angRepr = psi.GetRepresentation().GetRepresentation(angularRank)
	angSize = psi.GetData().shape[angularRank]

	for angIdx in range(angSize):
		curL = angRepr.Range.GetLmIndex(angIdx).l
		eigAngIdx = lmIdxList.index(LmIndex(curL, m))
		curE = eigenValues[eigAngIdx]
		curV = eigenVectors[eigAngIdx]

		boundIdx = where(curE < boundThreshold)[0]

		#Get projection on eigenstates
		psiSlice = psi.GetData()[angIdx, :]
		overlapPsi = dot(overlap, psiSlice)
		proj = dot(conj(curV[:,boundIdx].transpose()), overlapPsi)

		#Interpolate to get equispaced dP/dE
		boundDistr.append( proj )
		boundE.append( curE[boundIdx] )
		boundV.append(curV[:, boundIdx])
		boundTotal += sum(abs(proj)**2)

	return boundE, boundV, boundDistr, boundTotal


class EigenstateAnalysis:
	"""Perform analysis of a wavefunction in terms of eigenstates.
	"""
	
	def __init__(self, conf):
		self.Config = conf
		self.BoundThreshold = 0.0
		self.m = 0
		
	def Setup(self):
		#Setup pyprop problem
		self.Problem = pyprop.Problem(self.Config)
		self.Problem.SetupStep()
		
		#Calculate eigenstates
		E, V, angIdx, lmIdx = SetupRadialEigenstates(self.Problem)
		self.Eigenvalues = E
		self.EigenVectors = V
		self.AngularIndices = angIdx
		self.LMIndices = lmIdx
		
		#Setup overlap matrix
		self.Overlap =  SetupOverlapMatrix(self.Problem)
		
	
	def CalculateBoundProbability(self, psi):
		"""Calculate norm squared of bound part of psi
		"""
		dummy, dummy, dummy, boundTotal = \
			CalculateBoundDistribution(psi, self.Eigenvalues, self.EigenVectors, \
									self.LMIndices, self.m, self.Overlap, boundThreshold=self.BoundThreshold)
			
		return boundTotal
	
	
	def CalculateEnergyDistribution(self, psi):
		raise NotImplementedError("Not implemented yet!")
	