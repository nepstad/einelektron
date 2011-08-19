"""
EinPartikkel
============

The EinPartikkel package provides fundamental functionality related to one-
particle-type problems, such as potentials, analysis and eigenpair calculations.

"""
import sys
import os
import os.path as path

#Add Pyprop location to path
#PypropLocation = os.path.realpath("%s/../pyprop" % \
#		os.path.split(os.path.realpath(__file__))[0])
#sys.path.append(PypropLocation)

#Contains function and class name from package that might appear in config
#files. To be passed to Pyprop for resolving during loading of said files.
#from utils import ProjectNamespace

__all__ = ["analysis", "core", "eigenvalues", "utils", "quantumnumbers",\
"namegenerator"]


#Find out whether we running from a source or installed einpartikkel
IsRunningFromSource = True #not path.exists(path.join(__path__[0], "core", "libcore.so"))

if IsRunningFromSource:
	#find out which build kind we are
	BuildKind = os.environ.get("PYPROP_BUILD_KIND", "default")
	BuildPath = path.abspath(path.join(__path__[0], "..", "build", BuildKind, "einpartikkel"))


