
"""
eigenstates
===========

Basic eigenstate operations are provided. These are:

	1)Load/store eigenstates from/to disk.
	2)Iterate over bound/continuum states.


"""
from __future__ import with_statement
import os, errno
import tables
from numpy import array, where
import pyprop
#import ..
from ..eigenvalues.eigenvalues import SetupRadialEigenstates

#WTF is this?  
#@RegisterAll
class Eigenstates(object):
    def __init__(self, conf):
	"""
	Constructor method. Based on a config object, the eigenstates and eigenvalues
	are calculated and saved. 
	"""
	#Config object.
	self.Config = conf

	#Generate file name.
	self.FileName = self.NameGen(conf)
	
	#Store states if they have not been calculated/saved yet.
	if os.path.exists(self.FileName):
	    #Load eigenstates et al.
	    self.LoadEigenstates()
	else:
	    #Create pyprop object.
	    prop = pyprop.Problem(conf)
	    prop.SetupStep()
	    
	    #Retrieve lmax.
	    lmax = conf.AngularRepresentation.index_iterator.lmax
	    
	    #Calculate eigenstates.
	    E, V, angIdx, lmIdx = SetupRadialEigenstates(prop, mList = range(-lmax,lmax+1))
	    
	    #Saving states.
	    self.SaveEigenstates(E, V, angIdx, lmIdx)
	    
	    #Load eigenstates et al:
	    self.LoadEigenstates()

    
    def GetEigenstate(self, quantumNumber):
		"""
		Yet to be implemented.
		"""
		raise NotImplementedError("Not implemented yet!")
    
    def SaveEigenstates(self, eigenValues, eigenVectors, angIdxList, lmIdxList):
	"""
	SaveEigenstates(self, eigenValues, eigenVectors, angIdxList, lmIdxList)

	Saves eigenstates/values and corresponding lists of indices and quantum numbers.
	The format is meant o be read by LoadEigenstates() in this class. The data is stored in 
	a binary HDF5 format. The file name is an attribute of the class.

	Parametres
	----------
	eigenValues : list of 1D arrays of eigenvalues.
	eigenVectors : list of 2D arrays of eigenstates.
	angIdxList : integer list, showing which angular indices have been included.
	lmIdxList : list of {l,m} objects, showing wht quantum number each angular index corresponds to.

	"""
	#Check that folder(s) exists, make them if not
	filePath = os.path.dirname(self.FileName)
	if not os.path.exists(filePath):
		os.makedirs(filePath)
		try:
			os.makedirs(filePath)
		except OSError as exc: # Python >2.5
			if exc.errno == errno.EEXIST:
				pass
			else:
				raise

	f = tables.openFile(self.FileName, "w")
	try:
	    #Saving the lists of indices and results.
	    f.createArray("/", "Eigenvalues", eigenValues) 
	    f.createArray("/", "Eigenstates", eigenVectors)
	    f.createArray("/", "AngularIndex", angIdxList)
	    
	    #Saving config object.
	    f.setNodeAttr("/","Config",self.Config.cfgObj)

	    #Saving lm index list.
	    #Create array from object list. (Must be reversed in LoadEig...)
	    lm_list = []
	    for lm in lmIdxList:
		lm_list.append([lm.l, lm.m])
	    
	    f.createArray("/", "lm_list", lm_list)

	finally:
	    f.close()
	
    def LoadEigenstates(self):
	"""
	LoadEigenstates()

	Loads eigen-information that was saved by SaveEigenstates().
	"""
	f = tables.openFile(self.FileName, "r")
	try:
	    #Saving the lists of indices and results.
	    self.EigenValues = f.root.Eigenvalues[:]
	    self.EigenVectors = f.root.Eigenstates[:]
	    self.AngularIndices = f.root.AngularIndex[:]
	    
	    #Saving lm index list.
	    LMIndices = []
	    lmIdx = f.root.lm_list[:] 
	    for lm in lmIdx:
		LMIndices.append(pyprop.core.LmIndex(lm[0], lm[1]))
	    self.LMIndices = LMIndices
	finally:
	    f.close()

    def IterateStates(self, threshold):
	"""
	IterateStates(self, threshold)
	
	Iterate over all states with energies over threshold.

	Parametres
	----------
	threshold : float, the energy over which you want to look at the eigenstats.

	Example
	-------
	>>> conf = pyprop.Load("config.ini")
	>>> threshold = ionisation_probability = 0
	>>> my_eigenstates = Eigenstates(conf)
	>>> for angIdx, E, V, l, m in my_eigenstates.IterateStates(threshold):
	>>>	ionisation_probability += projection_method(my_wavefunction, V)	 
	"""
	for angIdx in self.AngularIndices:
	    curE = array(self.EigenValues[angIdx])
	    curV = array(self.EigenVectors[angIdx])
	    l = self.LMIndices[angIdx].l
	    m = self.LMIndices[angIdx].m
	    
	    #Filter out unwanted energies.
	    idx = where(curE > threshold)
	    curE = curE[idx]
	    curV = curV[:,idx]
	    yield angIdx, curE, curV, l, m

    def IterateBoundStates(self, threshold):
	"""
	IterateBoundStates(threshold)

	Iterate over all states with energies below threshold.	
	
	Parametres
	----------
	threshold : float, the energy UNDER which you want to look at the eigenstats.

	Example
	-------
	>>> conf = pyprop.Load("config.ini")
	>>> threshold = 0
	>>> my_eigenstates = Eigenstates(conf)
	>>> for angIdx, E, V, l, m in my_eigenstates.IterateStates(threshold):
	>>>	print angIdx, l, m
	 
	    0 0 0
	    1 1 -1
	    2 1 0
	    3 1 1
	    4 2 -2
	    . . .
	    . . .
	    . . .

	"""
	for angIdx in self.AngularIndices:
	    curE = array(self.EigenValues[angIdx])
	    curV = array(self.EigenVectors[angIdx])
	    l = self.LMIndices[angIdx].l
	    m = self.LMIndices[angIdx].m
	    
	    #Filter out unwanted energies.
	    idx = where(curE < threshold)
	    curE = curE[idx]
	    curV = curV[:,idx]
	    yield angIdx, curE, curV, l, m



    def NameGen(self, conf):
	"""
	fileName = NameGen(self, conf)

	Returns a tailored file name.

	Parametres
	----------
	conf : config object.

	Returns
	-------
	fileName : string, excellent file name.

	Notes
	-----
	This should be done better. This is a temporary arrangement. (I hope.)
	Assumes a folder named Eigenstates exists...
	"""
	#Folder
	my_name = "Eigenstates/"
	
	#Grid characteristics
	radialPostfix = "_".join(GetRadialPostfix(conf))
	angularPostfix = "_".join(GetAngularPostfix(conf))
	
	my_name += conf.Propagation.grid_potential_list[2] + "_"
	my_name += radialPostfix + angularPostfix
	
	return my_name




def GetRadialPostfix(conf):
    """
    GetRadialPostfix(conf)

    Returns a "unique" list of strings string identifying the radial grid
    implied by the specified args
    
    Parametres
    ----------
    conf : config object.
    """
    cfg = conf.RadialRepresentation

    gridType = cfg.bpstype
    postfix = ["grid", gridType, "xmax%i" % cfg.xmax, "xsize%i" % cfg.xsize, "order%i" % cfg.order]
    if gridType == "linear":
	    pass
    elif gridType == "exponentiallinear":
	    postfix.append("xpartition%i" % cfg.xpartition)
	    postfix.append("gamma%.1f" % cfg.gamma)
    elif gridType == "exponential":
	    postfix.append("gamma%.1f" % cfg.gamma)

    return postfix


def GetAngularPostfix(conf):
    """
    GetAngularPostfix(conf)

    Returns a "unique" list of strings string identifying the angular grid
    implied by the specified args
    
    Parametres
    ----------
    conf : config object.
    """
    postfix = ["angular"]
    postfix += ["lmax%i" % conf.AngularRepresentation.index_iterator.lmax]
    return postfix

