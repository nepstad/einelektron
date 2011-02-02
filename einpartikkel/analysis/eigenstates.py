
"""
eigenstates
===========

Basic eigenstate operations are provided. These are:

	1)Load/store eigenstates from disk


"""
from __future__ import with_statement
import tables
from numpy import array, where
import pyprop
import os.path
#import ..
from ..eigenvalues.eigenvalues import SetupRadialEigenstates

#WTF is this?  
#@RegisterAll
class Eigenstates(object):
    def __init__(self, conf):
	"""
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
	    
	    #Load eigenstates et al.
	    self.LoadEigenstates()

    
    def GetEigenstate(self, quantumNumber):
	raise NotImplementedError("Not implemented yet!")
    
    def SaveEigenstates(self, eigenValues, eigenVectors, angIdxList, lmIdxList):
	"""

	"""
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

	Iterate over all states with energies over threshold.
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

	Iterate over all states with energies below threshold.
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




#@RegisterAll
def GetRadialPostfix(conf):
    """
    Returns a "unique" list of strings string identifying the radial grid
    implied by the specified args
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


#@RegisterAll
def GetAngularPostfix(conf):
    """
    Returns a "unique" list of strings string identifying the angular grid
    implied by the specified args
    """
    postfix = ["angular"]
    postfix += ["lmax%i" % conf.AngularRepresentation.index_iterator.lmax]
    return postfix

