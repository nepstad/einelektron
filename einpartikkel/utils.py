"""
Utils
=====

Misc utilities are put here for now.
"""

import sys
import logging
import inspect
import pyprop.core

ProjectNamespace = []


def RegisterAll(obj):
	"""
	Function decorator for registering 'obj' in containing module's
	__all__ variable. 'obj' could be a class, or a function, as long
	as it has a __name__ and __module__ attribute.

	Based on ideas from:
	  http://groups.google.com/group/comp.lang.python/msg/11cbb03e09611b8a
	  http://code.activestate.com/recipes/576993/
	"""
	
	all = sys.modules[obj.__module__].__dict__.setdefault('__all__', [])
	if obj.__name__ not in all:  # Prevent duplicates if run from an IDE.
		all += [obj.__name__]
	return obj
RegisterAll(RegisterAll)


@RegisterAll
def RegisterProjectNamespace(obj):
	"""
	Function decorator for registering 'obj' in the ProjectNamespace
	variable for EinElektron, to be used by Pyprop for resolving config files 
	content.

	Based on ideas from:
	  http://groups.google.com/group/comp.lang.python/msg/11cbb03e09611b8a
	  http://code.activestate.com/recipes/576993/
	"""
	global ProjectNamespace
	# Prevent duplicates if run from an IDE.
	if obj.__name__ not in ProjectNamespace:  		
		ProjectNamespace += [obj]
	return obj


@RegisterAll
def UpdatePypropProjectNamespace(pypropProjNamespace):
	"""
	Update Pyprop project namespace with objects from this (EinPartikkel) package

	Input
	-----
	pypropProjNamespace: Pyprop project namespace

	Output
	------
	None
	"""
	for ref in ProjectNamespace:
		if ref not in pypropProjNamespace.keys():
			pypropProjNamespace[ref.__name__] = ref


