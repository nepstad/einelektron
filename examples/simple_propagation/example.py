from numpy import sin, cos, pi
import sys
sys.path.append("../..")
import einpartikkel
import pyprop

from einpartikkel.propagation.propagate import Propagate
from einpartikkel.propagation.tasks import ComputeAtomicInitialState, ProgressReport, DisplayGMRESError
from einpartikkel import quantumnumbers

import einpartikkel.core.indexiterators
import einpartikkel.core.preconditioner
from einpartikkel.utils import UpdatePypropProjectNamespace
UpdatePypropProjectNamespace(pyprop.ProjectNamespace)


def RunPropagation():
	#Load config
	configFile = "config.ini"
	conf = pyprop.Load(configFile)
	
	#Setup propagation tasks
	qnum = quantumnumbers.HydrogenicQuantumNumbers(1,0,0)
	tasks = [ComputeAtomicInitialState(qnum), ProgressReport(), DisplayGMRESError()]

	#Setup propagation
	prop = Propagate(conf, tasks, 50)

	#Run propagation
	prop.run()

	return prop
	


def LaserFunctionSimpleLength(conf, t):
	if 0 <= t < conf.pulse_duration:
		curField = conf.amplitude;
		curField *= sin(t * pi / conf.pulse_duration)**2;
		curField *= cos(t * conf.frequency);
	else:
		curField = 0
	return curField
pyprop.ProjectNamespace["LaserFunctionSimpleLength"] = LaserFunctionSimpleLength

