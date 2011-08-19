"""
Core
====

Core functionality for EinPartikkel calculations, such potentials.

-Potential implementations (c++)
-Time functions for time-dependent potentials (python)
-Preconditioner setup for one-particle (3D) problems

"""
__all__ = ["analysis", "namegenerator", "eigenstates", "coulombwaves"]

import sys
import os

import pyprop
if pyprop.IsRunningFromSource:
    sys.path.append(os.path.join(pyprop.BuildPath, "../modules", "einpartikkel", "core"))
import libeinelektroncore

from libeinelektroncore import *
