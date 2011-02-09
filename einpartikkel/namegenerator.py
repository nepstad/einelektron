"""
namegenerator
-------------

Tools to autogenerate names from config objects

"""
from utils import RegisterAll

@RegisterAll
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


@RegisterAll
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

