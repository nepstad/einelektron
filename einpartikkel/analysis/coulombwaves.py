"""
coulombwaves
============

Basic Coulomb wave operations are provided. These are:

	1) construct Coulomb waves
	2) store to disk
	3) load from disk
	4) nice class representation of Coulomb waves

"""
import sys
import os, errno
import tables
from numpy import array, where, r_, sqrt, zeros, array, abs, dot, conj, pi, \
	double, complex, iterable, diff
import pyprop
from einpartikkel.eigenvalues.eigenvalues import SetupRadialEigenstates
from einpartikkel.utils import RegisterAll
from einpartikkel.namegenerator import GetRadialPostfix, GetAngularPostfix
from above import SetRadialCoulombWave

@RegisterAll
class CoulombWaves(object):
	"""
	"""
	_DataItemAttrName = "DataItems"
	_MetaDataItemAttrName = "MetaDataItems"
	
	def __init__(self, conf, Z, Emax, dE):
		"""Return a Coulomb wave object
		"""
		#Config object.
		self.Config = conf

		self.Z = Z
		self.MaxEnergy = Emax
		self.EnergyResolution = dE

		#Data structure for Coulomb waves etc
		self._Data = {}
		self._MetaData = {}

		md = self._MetaData
		md["Z"] = self.Z
		md["Emax"] = self.MaxEnergy
		md["dE"] = self.EnergyResolution

		#Generate file name.
		self.FolderName = "CoulombWaves" #TODO: this should be generated somehow
		self.FileName = CoulombWaveNameGenerator(self.Config, self.FolderName,\
				self.Z, dE, Emax)

		#Get logger handle
		self.Logger = pyprop.GetClassLogger(self)

		self._SetupComplete = False
	

	def SetData(self, data):
		self._Data = data
		self._SetupComplete = True


	@classmethod
	def FromFile(cls, filename):
		"""Setup CoulombWaves object from serialized (file) data

		Takes a filename as argument.
		"""
		if iterable(filename):
			obj = cls.FromMultipleFiles(cls, filename)
		else:
			obj = cls.FromMultipleFiles(cls, [filename])

		return obj


	@classmethod
	def FromMultipleFiles(cls, fileList):
		"""Setup CoulombWaves object from serialized (file) data

		Takes a list of file names as argument
		"""
		logger = pyprop.GetFunctionLogger()
		logger.info("Setting up Coulomb waves from serialized data...")

		d = {}
		metaDataItems = {}
		for filename in fileList:
			with tables.openFile(filename, "r") as f:
				#Load data
				dataItemNames = f.getNodeAttr("/", cls._DataItemAttrName)
				for dataItem in dataItemNames:
					d[dataItem] = f.getNode("/", dataItem)[:]

				#Load metadata
				metaDataItemNames = f.getNodeAttr("/", cls._MetaDataItemAttrName)
				for mdItem in metaDataItemNames:
					metaDataItems[mdItem] = f.getNodeAttr("/", mdItem)
				
			#Load config object.
			config = pyprop.LoadConfigFromFile(filename, "/")
			

		#Instantitate CoulombWave object and set data
		obj = cls(config, **metaDataItems)
		obj.SetData(d)
		
		logger.info("Loading complete.")

		return obj


	def Setup(self):
		"""Calculate the Coulomb waves specified by config object
		
		"""
		
		#Check if we have already calculated coulomb waves
		if self._SetupComplete:
			return

		d = self._Data

		#Retrieve list of angular momenta (l)
		d["AngularMomenta"] = [it.l for it in self.Config.AngularRepresentation.index_iterator]

		#Create wavefunction
		self.Psi = pyprop.CreateWavefunction(self.Config)
		
		#Calculate eigenstates
		for l in d["AngularMomenta"]:
			E, cw = SetupRadialCoulombStatesEnergyNormalized(self.Psi, self.Z, \
				l, self.MaxEnergy, self.EnergyResolution)

			d["Energies_l%i" % l] = E
			d["CoulombWaves_l%i" % l] = cw

		self._SetupComplete = True
    
    
	def Serialize(self):
		"""Store Coulomb waves on disk

		Store Coulomb waves in a HDF5 file with autogenerated name

		"""
		assert (self._SetupComplete)

		#Check that folder(s) exists, make them if not
		filePath = os.path.dirname(self.FileName)
		if not os.path.exists(filePath):
			os.makedirs(filePath)
			try:
				os.makedirs(filePath)
			except OSError as exc:
				if exc.errno == errno.EEXIST:
					pass
				else:
					raise

		self.Logger.info("Now saving Coulomb waves to disk...")
		with tables.openFile(self.FileName, "w") as f:
			#Store Coulomb wave data
			for dataItemName, dataItemValue in self._Data.iteritems():
				f.createArray("/", "%s" % dataItemName, dataItemValue)

			#Store metadata
			for dataItemName, dataItemValue in self._MetaData.iteritems():
				f.setNodeAttr("/", dataItemName, dataItemValue)
			
			#Store data item names
			f.setNodeAttr("/", self._DataItemAttrName, self._Data.keys())
			f.setNodeAttr("/", self._MetaDataItemAttrName, self._MetaData.keys())

		#Saving config object.
		pyprop.serialization.SaveConfigObject(self.FileName, "/", self.Config)

		self.Logger.info("Save complete.")


	def IterateStates(self, threshold):
		"""
		IterateStates(self, threshold)
		
		Iterate over all states with energies over threshold.

		Parametres
		----------
		threshold : (float) lower energy cutoff

		"""
		assert (self._SetupComplete)
		#enList = self._Data["Energies"]
		#cwList = self._Data["CoulombWaves"]
		#for curL, (E, curCW) in enumerate(zip(enList, cwList)):
		idxIt = self.Config.AngularRepresentation.index_iterator
		for angIdx, lmIdx in enumerate(idxIt):
			#Get current l and m
			l = lmIdx.l
			m = lmIdx.m
			
			#Get energies and Coulomb waves
			curE = self._Data["Energies_l%i" % l]
			curCW = self._Data["CoulombWaves_l%i" % l]
			
			#Get energies and Coulomb waves
			#curE = array(enList[l])
			#curCW = array(cwList[l])

			#Filter out unwanted energies.
			idx = where(curE > threshold)
			filteredE = curE[idx]
			filteredCW = curCW[:,idx]
			#filteredCW = curCW[idx,:]
			
			#Test:normalize
			k = sqrt(2*filteredE)
			dE = diff(filteredE)[0]
			factor = sqrt(2*dE/(pi*k))


			yield angIdx, filteredE, filteredCW * factor, l, m


@RegisterAll
def CoulombWaveNameGenerator(conf, folderName, Z, dE, Emax):
	"""fileName = CoulombWaveNameGenerator(conf, folderName)

	Returns a generated file name.

	Input
	-----
	conf: Pyprop config
	folderName: (string) name of folder where Coulomb waves are stored
	Z: (int) Coulomb charge
	dE: (double) energy spacing of Coulomb waves
	Emax: (double) highest Coulomb wave energy

	Returns
	-------
	fileName : (string) generated file name.

	"""

	def getAngularPostfix(conf):
		lList = [it.l for it in conf.AngularRepresentation.index_iterator]
		return ["l%s" % lList]

	#Coulomb wave characteristics
	radialPostfix = "_".join(GetRadialPostfix(conf))
	#angularPostfix = "_".join(GetAngularPostfix(conf))
	angularPostfix = "_".join(getAngularPostfix(conf))
	chargePostfix = "Z%s" % Z
	energyPostfix = "dE%s_Emax%s" % (dE, Emax)
	
	filename = folderName + "/" 
	filename += "_".join([radialPostfix, angularPostfix, \
		chargePostfix, energyPostfix])
	filename += ".h5"
	
	return filename


#---------------------------------------------------------------------------------------
#            Coulomb Wave Analysis
#---------------------------------------------------------------------------------------
def SetupRadialCoulombStatesEnergyNormalized(psi, Z, l, Emax, dE):
	"""Calculate radial Coulomb states with energy normalization

	Input
	-----
	psi:  (wavefunction) a wavefunction with the desired representation
	Z:    (int) the Coulomb charge
	l:    (int) angular momentum
	Emax: (double) desired max energy of Coulomb waves
	dE:   (double) energy resolution

	"""

	logger = pyprop.GetFunctionLogger()
	E = r_[dE:Emax:dE]
	k = sqrt(E*2)
	
	bspline = psi.GetRepresentation().GetRepresentation(1).GetBSplineObject()
	#l = array(psi.GetRepresentation().GetGlobalGrid(0), dtype=int)

	#Setup Radial Waves for given l
	logger.info("Generating Coulomb waves for Z = %s, l = %s..." % (Z, l))
	f = GetRadialCoulombWaveBSplines
	def cwFunc(curK):
		idx = curK[0]
		k = curK[1]
		#Normalization factor
		factor = sqrt(2*dE/(pi*k))
		return factor * f(Z, int(l), k, bspline)
	result = array(map(cwFunc, enumerate(k)))

	coulWaves = result.transpose()
	energies = E

	return energies, coulWaves


def GetRadialCoulombWaveBSplines(Z, l, k, bsplineObj):
	"""Calculate radial Coulomb wave expanded in B-splines

	Input
	-----
	Z: (int) Coulomb charge
	l: (int) angular momentum
	k: (double) Coulomb wave momentum
	bsplineObj: BSpline instance

	"""

	#Setup buffers for the calculation
	def setupBuf():
		r = bsplineObj.GetQuadratureGridGlobal()
		wav = zeros(len(r), dtype=double)
		coeff = zeros(bsplineObj.NumberOfBSplines, dtype=complex)
		return r, wav, coeff
	r, wavBuf, coeff = setupBuf()

	#Get the Coulomb function in grid space
	def setCoulombwave():
		SetRadialCoulombWave(Z, l, k, r, wavBuf)
	setCoulombwave()

	#Copy coulomb waves to a complex buffer
	def copyBuf():
		cplxWav = array(wavBuf, dtype=complex)
		return cplxWav
	cplxWavBuf = copyBuf()

	#get bspline coeffs
	def expandFunc():
		bsplineObj.ExpandFunctionInBSplines(cplxWavBuf, coeff)
	expandFunc()

	return coeff

