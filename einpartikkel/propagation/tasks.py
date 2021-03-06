"""
PropagationTasks
================

Defines standard tasks to be performed during propagation.

A list of instantiated propagation tasks should be supplied to the Propagate
class. 

Custom propagation tasks should derive from the PropagationTask baseclass,
defined here.

"""

from __future__ import with_statement
import os.path
import time
import tables
import pypar
import pyprop
from pyprop import PrintOut
from pyprop.pyproplogging import GetClassLogger, GetFunctionLogger
from ..utils import RegisterAll
from ..eigenvalues import eigenvalues


def CreatePath(absFileName):
	"""Create directories in abspath
	
	"""
	logger = GetFunctionLogger()
	if pyprop.ProcId == 0:
		filePath = os.path.dirname(absFileName)
		if not os.path.exists(filePath) and len(filePath) > 0:
			logger.debug("Creating folder: %s" % filePath)
			os.makedirs(filePath)
	pypar.barrier()
	

class PropagationTask:
	def __init__(self):
		raise NotImplementedError("Please implement in derived class")
	
	def setupTask(self, prop):
		raise NotImplementedError("Please implement in derived class")


	def callback(self, prop):
		raise NotImplementedError("Please implement in derived class")

	def postProcess(self, prop):
		raise NotImplementedError("Please implement in derived class")


@RegisterAll
class ProgressReport(PropagationTask):
	"""
	Print some progress information during propagation, and store it
	in a HDF5-file when finished (time, norm and projection on initial state)
	"""

	def __init__(self, storeProgressInfo=True):
		self.StoreProgressInfo = storeProgressInfo
		self.Logger = GetClassLogger(self)
		self.StartTime = -1
		self.InitialPsi = None
		self.ProgressItems = {}

	def setupTask(self, prop):
		self.Logger.info("Setting up task...")
		self.StartTime = time.time()
		self.InitialPsi = prop.psi.Copy()
	
		#check if output dir exist, create if not
		if self.StoreProgressInfo:
			self.OutputFileName = prop.Config.Names.output_file_name
			CreatePath(self.OutputFileName)
		
		#Report items
		self.ProgressItems["SampleTimes"] = []
		self.ProgressItems["Norm"] = []
		self.ProgressItems["InitialCorrelation"] = []

	def callback(self, prop):
		t = prop.PropagatedTime
		T = prop.Duration + prop.StartTime
		norm = prop.psi.GetNorm()
		corr = abs(prop.psi.InnerProduct(self.InitialPsi))**2
		eta = self.__EstimateETA(prop)
		self.ProgressItems["SampleTimes"] += [t]
		self.ProgressItems["Norm"] += [norm]
		self.ProgressItems["InitialCorrelation"] += [corr]
		#FormatDuration = lambda t: time.strftime("%Hh %Mm %Ss", time.gmtime(t))
		#FormatDuration = lambda t: time.strftime("%dd %Hh %Mm %Ss", time.gmtime(t))
		PrintOut("t = %.2f / %.2f; N = %.15f; Corr = %.12f, ETA = %s" % (t, T, norm, corr, self._FormatDuration(eta)))

	def postProcess(self, prop):
		"""
		Store problem information collected during propagation
		"""
		if self.StoreProgressInfo and (pyprop.ProcId == 0):
			with tables.openFile(self.OutputFileName, "a") as h5file:
				for itemName, itemVal in self.ProgressItems.iteritems():
					if itemName in h5file.root:
						h5file.removeNode(h5file.root, itemName, recursive=True)
					h5file.createArray("/", itemName, itemVal)


	def __EstimateETA(self, prop):
		"""
		Estimates remaining time before propagation is finished
		"""
		curTime = time.time() - self.StartTime
		#totalTime = (curTime / prop.PropagatedTime) * prop.Duration
		totalTime = (curTime / (prop.PropagatedTime - prop.StartTime)) * prop.Duration
		eta = totalTime - curTime
		return eta

	def _FormatDuration(self, t):
		days, remainder = divmod(t, 24 * 60 * 60)
		hours, remainder = divmod(remainder, 60 * 60)
		minutes, seconds = divmod(remainder, 60)
		if days > 0:
			timeStr = "%id %ih" % (days, hours)
		elif hours > 0:
			timeStr = "%ih %im" % (hours, minutes)
		else:
			timeStr = "%im %is" % (minutes, seconds)
			
		return timeStr


class DisplayGMRESError(PropagationTask):
	"""
	Print GMRES solver error at each callback
	"""
	def __init__(self):
		pass
	
	def setupTask(self, prop):
		pass
	
	def callback(self, prop):
		PrintOut(prop.Propagator.Solver.GetErrorEstimateList())
		
	def postProcess(self, prop):
		pass


class SaveWavefunction(PropagationTask):
	"""
	Save wavefunction after propagation, and, if specified, for each
	callback during propagation
	"""

	def __init__(self, storeDuringPropagation):
		self.StoreDuringPropagation = storeDuringPropagation
		self.Counter = 0

	def setupTask(self, prop):
		#get output filename
		self.OutputFileName = prop.Config.Names.output_file_name

		#check if output dir exist, create if not
		CreatePath(self.OutputFileName)

		#store the initial wavefunction
		prop.SaveWavefunctionHDF(self.OutputFileName, "/initialWavefunction")

	def callback(self, prop):
		if self.StoreDuringPropagation:
			#create unique filename
			filename = "%s_%03i.h5" % (self.OutputFileName.strip(".h5"), self.Counter)
			
			#store current wavefunction and propagation time
			prop.SaveWavefunctionHDF(filename, "/wavefunction")
			if pyprop.ProcId == 0:
				with tables.openFile(filename, "r+") as h5:
					h5.setNodeAttr("/wavefunction", "prop_time", prop.PropagatedTime)
			pypar.barrier()

			self.Counter += 1

	def postProcess(self, prop):
		#store the final wavefunction
		prop.SaveWavefunctionHDF(self.OutputFileName, "/wavefunction")


class ComputeAtomicInitialState(PropagationTask):
	"""
	Diagonalize problem hamiltonian to determine eigenstates, and then
	set initial wavefunction to one of these eigenstates.

	"""

	def __init__(self, initialStateQuantumNumbers):
		self.QuantumNumbers = initialStateQuantumNumbers

	def setupTask(self, prop):
		"""Calculate bound state and set prop.psi equal specified one.
		"""
		E, V, angIdxList, lmIdxList = eigenvalues.SetupRadialEigenstates(prop, potentialIndices=[0], mList=[self.QuantumNumbers.m])
		eigenvalues.SetRadialEigenstate(prop.psi, V, angIdxList, self.QuantumNumbers)

	def callback(self, prop):
		pass

	def postProcess(self, prop):
		pass


