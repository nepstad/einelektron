from numpy import sin, cos, pi
import sys
sys.path.append("../..")
import pyprop
from pyprop.modules import einpartikkel

#Imports needed by config
from pyprop.core.representation import CombinedRepresentation_2
from pyprop.modules.discretization.bsplines import BS

from einpartikkel.propagation.propagate import Propagate
from einpartikkel.propagation.tasks import ComputeAtomicInitialState, ProgressReport, DisplayGMRESError, \
		SaveWavefunction
from einpartikkel import quantumnumbers

import einpartikkel.core.indexiterators
import einpartikkel.core.preconditioner

from libeinpartikkelcore import *


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



#Laser function velocity gauge
def LaserFunctionSimpleVelocity(conf, t):
	"""
		Laser function velocity gauge.
	"""
	if 0 <= t < conf.pulse_duration:
		curField = conf.amplitude / conf.frequency;
		curField *= sin(t * pi / conf.pulse_duration)**2;
		curField *= sin(t * conf.frequency + conf.phase);
	else:
		curField = 0
	return curField


#Laser function length gauge
def LaserFunctionSimpleLength(conf, t):
	"""
		Laser function length gauge.
	"""
	if 0 <= t < conf.pulse_duration:
		curField1 = pi / (conf.pulse_duration * conf.frequency);
		curField1 *= sin(2 * pi * t / conf.pulse_duration);
		curField1 *= sin(conf.frequency * t + conf.phase);

		curField2 = sin(pi * t / conf.pulse_duration);
		curField2 *= sin(pi * t / conf.pulse_duration);
		curField2 *= cos(conf.frequency * t + conf.phase);

		curField = -conf.amplitude * (curField1 + curField2);

	else:
		curField = 0
	return curField



#Laser functions in each direction (Velocity gauge)
def LaserFunctionSimpleVelocity_X(conf, t):
	return LaserFunctionSimpleVelocity(conf, t)

def LaserFunctionSimpleVelocity_Y(conf, t):
	return LaserFunctionSimpleVelocity(conf, t)

def LaserFunctionSimpleVelocity_Z(conf, t):
	return LaserFunctionSimpleVelocity(conf, t)


#Laser function in each direction (Length gauge)
def LaserFunctionSimpleLength_X(conf, t):
	return LaserFunctionSimpleLength(conf, t)

def LaserFunctionSimpleLength_Y(conf, t):
	return LaserFunctionSimpleLength(conf, t)

def LaserFunctionSimpleLength_Z(conf, t):
	return LaserFunctionSimpleLength(conf, t)

