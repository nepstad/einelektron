import pypar
import sys
sys.path.append("../..")
import einpartikkel
import pyprop

from einpartikkel.propagation.propagate import Propagate
from einpartikkel.propagation.tasks import ComputeAtomicInitialState, ProgressReport, DisplayGMRESError, \
		SaveWavefunction, PropagationTask
from einpartikkel import quantumnumbers
from einpartikkel.eigenvalues.eigenvalues_iter import \
	SaveEigenpairsSolverShiftInvert, FindEigenvaluesInverseIterationsPiram

import einpartikkel.core.indexiterators
import einpartikkel.core.preconditioner
from einpartikkel.utils import UpdatePypropProjectNamespace
UpdatePypropProjectNamespace(pyprop.ProjectNamespace)


def RunPropagation():
	"""Perform propagation defined by 'config.ini' file.

	"""

	#Load config
	configFile = "config.ini"
	conf = pyprop.Load(configFile)
	
	#Setup propagation tasks. Initial state = ground state
	qnum = quantumnumbers.HydrogenicQuantumNumbers(2,1,0)
	tasks = [ComputeAtomicInitialState(qnum), ProgressReport(), DisplayGMRESError(), SaveWavefunction(False)]

	#Setup propagation
	prop = Propagate(conf, tasks, 100)
	prop.preprocess()

	#Run propagation
	prop.run()

	return prop
	
	
def FindGroundstate(configFile="groundstate.ini"):
	"""Find groundstate using imaginary time propagation
	
	"""
	timers = pyprop.Timers()

	#Load config
	timers["Setup"].Start()
	conf = pyprop.Load(configFile)

	#Update config
	conf.SetValue("Names", "output_file_name", "./groundstate_fd.h5")
	conf.SetValue("Propagation", "timestep", -0.1j)
	conf.SetValue("Propagation", "duration", 50.0)
	conf.SetValue("Propagation", "renormalization", True)

	try:
		conf.Propagation.grid_potential_list.remove("Absorber")
	except:
		pass

	#Tasks
	tasks = [ProgressReport(), EnergyReport(), DisplayGMRESError()]#, SaveWavefunction(False)]

	#Setup propagation
	prop = Propagate(conf, tasks, 10)
	prop.preprocess()
	timers["Setup"].Stop()


	#Run propagation
	timers["Propagate"].Start()
	prop.run()
	timers["Propagate"].Stop()

	if pyprop.IsMaster():
		print
		print timers

	return prop

	
def FindEigenstates(configFile="groundstate.ini"):
	"""Find eigenstates
	
	"""
	#Load config
	conf = pyprop.Load(configFile)

	#Setup solver
	solver, siSolver, E = FindEigenvaluesInverseIterationsPiram(conf)

	#Find eigs
	SaveEigenpairsSolverShiftInvert(solver, siSolver)

	if pyprop.IsMaster:
		print E.real


class EnergyReport(PropagationTask):
	def __init__(self):
		pass

	def setupTask(self, prop):
		pass

	def preprocess(self, prop):
		pass

	def callback(self, prop):
		curEnergy =  prop.GetEnergyExpectationValue()
		pyprop.PrintOut("Energy = %s" % curEnergy)

	def postProcess(self, prop):
		pass


if __name__ == "__main__":
	#FindGroundstate()
	FindEigenstates()
	pypar.finalize()
