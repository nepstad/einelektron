import pyprop
from pyprop.debug import PrintMemoryUsage
import pyprop.modules.solvers.trilinos as trilinos
from ..utils import RegisterAll, RegisterProjectNamespace

@RegisterAll
class RadialPreconditioner:
	"""
	Preconditioner for GMRES for solving systems of the type

	(1)     (a H + b S) x = y  ->  x

	where H is the Hamiltonian, S is the overlap matrix
	and a and b are complex scaling factors (scalingH and scalingS).

	The preconditioner is given a set of tensor potentials which are
	diagonal in the angular rank (also called radial potentials).
	These potentials should approximate the Hamiltonian as well as
	possible, as the system (1) is solved exactly for the given radial
	potentials
	"""

	def __init__(self, psi):
		self.Rank = psi.GetRank()
		self.psi = psi

	def ApplyConfigSection(self, conf):
		self.OverlapSection = conf.Config.GetSection("OverlapPotential")
		self.PotentialSections = [conf.Config.GetSection(s) for s in conf.potential_evaluation]

	def SetHamiltonianScaling(self, scalingH):
		self.HamiltonianScaling = scalingH

	def SetOverlapScaling(self, scalingS):
		self.OverlapScaling = scalingS

	def GetHamiltonianScaling(self):
		return self.HamiltonianScaling

	def GetOverlapScaling(self):
		return self.OverlapScaling

	def Setup(self, prop):
		"""
		Set up a tensor potential for overlap potential and all other potentials
		and add them together, assuming they have the same layout
		ending up with a potential containing S + scalingH * (P1 + P2 + ...)

		The radial part of this potential is then converted to compressed col storage
		and factorized.
		"""

		#Setup overlap potential
		PrintMemoryUsage("Before Preconditioner Generate Potential (Overlap)")
		tensorPotential = prop.BasePropagator.GeneratePotential(self.OverlapSection)
		tensorPotential.PotentialData *= self.GetOverlapScaling()

		#Add all potentials to solver
		scalingH = self.GetHamiltonianScaling()
		for conf in self.PotentialSections:
			PrintMemoryUsage("Before Preconditioner Generate Potential (%s)" % conf)
			#Setup potential in basis
			potential = prop.BasePropagator.GeneratePotential(conf)
			if not tensorPotential.CanConsolidate(potential):
				raise Exception("Cannot consolidate potential %s with overlap-potential" % (potential.Name))

			#Add potential
			potential.PotentialData *= scalingH
			tensorPotential.PotentialData += potential.PotentialData
			del potential
		PrintMemoryUsage("After Preconditioner Generate Potentials")


		#Setup solvers
		tensorPotential.SetupStep(0.0)
		self.SetupRadialSolvers(tensorPotential)

	def SetupRadialSolvers(self, tensorPotential):
		raise NotImplementedException("Please Override")

	def Solve(self, psi):
		raise NotImplementedException("Please Override")

@RegisterAll
@RegisterProjectNamespace
class RadialPreconditionerIfpack(RadialPreconditioner):
	"""
	RadialPreconditioner using Ifpack (ILU) to
	approximately factorize the radial blocks
	"""

	def __init__(self, psi):
		RadialPreconditioner.__init__(self, psi)

	def ApplyConfigSection(self, conf):
		RadialPreconditioner.ApplyConfigSection(self, conf)
		self.Cutoff = conf.cutoff

	def SetupRadialSolvers(self, tensorPotential):
		PrintMemoryUsage("Before Ifpack Setup")

		#Setup the ILU preconditioner for each radial rank
		radialSolvers = []
		matrixCount = tensorPotential.PotentialData.shape[0]
		assert(tensorPotential.PotentialData.shape[0] == self.psi.GetData().shape[0])
		for i in range(matrixCount):
			vector = self.psi.GetData()[i,:]
			matrix = tensorPotential.PotentialData[i, :]

			solver = pyprop.createinstance.CreateInstanceRank("trilinos.IfpackRadialPreconditioner", 1)
			basisPairs = tensorPotential.BasisPairs[1:]
			solver.Setup(vector, matrix, basisPairs, self.Cutoff)
			radialSolvers.append(solver)

		self.RadialSolvers = radialSolvers
		PrintMemoryUsage("After Ifpack Setup")

	def Solve(self, psi):
		data = psi.GetData()

		angularCount = data.shape[0]
		if angularCount != len(self.RadialSolvers):
			raise Exception("Invalid Angular Count")

		for angularIndex, solve in enumerate(self.RadialSolvers):
			solve.Solve(data[angularIndex, :])



