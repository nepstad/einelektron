"""
analysis
========

Subpackage for performing various analysis tasks for EinPartikkel problems.

"""

__all__ = ["analysis", "namegenerator", "eigenstates", "coulombwaves"]
import sys
import os

import pyprop
if pyprop.IsRunningFromSource:
    sys.path.append(os.path.join(pyprop.BuildPath, "../modules", "einpartikkel", "analysis"))

import libeinpartikkelanalysis
from libeinpartikkelanalysis import SetRadialCoulombWave

