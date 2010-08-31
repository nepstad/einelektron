"""
Core
====

Core functionality for EinPartikkel calculations, such potentials.

-Potential implementations (c++)
-Time functions for time-dependent potentials (python)
-Preconditioner setup for one-particle (3D) problems

"""

#put boost::python-wrapped potentials into project namespace
from ..utils import ProjectNamespace, RegisterProjectNamespace
import above
from above import *
for key in above.__dict__.iterkeys():
	if not key.startswith("__"):
		RegisterProjectNamespace(eval(key))

__all__ = ["above", "preconditioner", "indexiterators"]
