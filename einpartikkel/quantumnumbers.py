import pyprop

from utils import RegisterAll
from pyprop.modules.discretizations.sphericalbasis import LmIndex

@RegisterAll
class HydrogenicQuantumNumbers(object):
	"""
	A simple class to handle quantum numbers for hydrogenic systems,
	that is (n,l,m)

	"""

	def __init__(self, n, l, m):
		assert (n >= 0)
		assert (l < n)
		assert (abs(m) <= l)
		self.n = n
		self.l = l
		self.m = m

	def __eq__(self, other):
		return (self.n == other.n) and (self.l == other.l) and (self.m == other.m)

	def GetRadialIndex(self):
		return self.n - self.l - 1

	def GetLmIndex(self):
		return LmIndex(self.l, self.m)
