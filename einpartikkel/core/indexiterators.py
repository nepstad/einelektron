import numpy
from numpy import unique
#Import CoupledIndex from core, and add a pretty print to it
def LmIndexIter(self):
	yield self.l
	yield self.m

from pyprop.modules.discretizations.sphericalbasis import LmIndex
LmIndex.__str__ = lambda self: "l=%i, m=%i" % (self.l, self.m)
LmIndex.__repr__ = lambda self: "LmIndex(%i, %i)" % (self.l, self.m)
LmIndex.__iter__ = LmIndex

from ..utils import RegisterAll, RegisterProjectNamespace

@RegisterAll
@RegisterProjectNamespace
class DefaultLmIndexIterator:
	"""
	Creates an iterator for giving all possible values of l and m to
	the SphericalHarmonicBasisRepresentation

	lmax gives the maximum l value of l, while m takes on all values
	allowed by |m| <= l

	This generator is most likely used in a configuration file like

	[AngularRepresentation]
	type = core.SphericalHarmonicBasisRepresentation
	index_iterator = DefaultLmIndexIterator(lmax=4)

	"""
	def __init__(self, lmax):
		self.lmax = lmax

	def __iter__(self):
		#Iterate through all permutations
		for curl in range(self.lmax + 1):
			for curm in range(-curl, curl + 1):
				yield LmIndex(curl,curm)

	def __repr__(self):
		return "%s(%s)" % \
			(self.__class__.__name__, self.lmax)


@RegisterAll
@RegisterProjectNamespace
class FixedMLmIndexIterator:
	"""
	Creates an iterator for giving all possible values of l and a given list
	of m values to the SphericalHarmonicBasisRepresentation.

	lmax gives the maximum l value of l.

	m should be a list of values, with the constraint |m| <= l

	This generator is most likely used in a configuration file like

	[AngularRepresentation]
	type = core.SphericalHarmonicBasisRepresentation
	index_iterator = DefaultLmIndexIterator(lmax=4, m=[-1,0,1])

	"""
	def __init__(self, lmax, m):
		self.lmax = lmax
		self.m = m

		assert(numpy.iterable(m))
		assert(abs(max(m)) <= lmax)

	def __iter__(self):
		#Iterate through all permutations
		for curm in self.m:
			for curl in range(abs(curm), self.lmax + 1):
				yield LmIndex(curl,curm)

	def __repr__(self):
		return "%s(%s, %s)" % \
			(self.__class__.__name__, self.lmax, self.m)
