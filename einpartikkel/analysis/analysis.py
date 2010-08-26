from pyprop.core import LmIndex
from numpy import dot, conj, where

def CalculateBoundDistribution(psi, eigenValues, eigenVectors, lmIdxList, m, overlap, boundThreshold=0):
	"""
	Calculate norm**2 of projection onto all bound states. The given eigenvectors
	are assumed to cover all different angular momenta, but only m=0 (i.e. no 
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
