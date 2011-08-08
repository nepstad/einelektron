import os.path
import os

APPNAME='pyprop-einpartikkel'
VERSION='1.0'

srcdir = '.'
blddir = 'build'

def set_options(opt):
	gr = opt.parser.add_option_group("pyprop build options")
	opt.pyprop_build_options = gr
	gr.add_option("--pyste", action="store_true", dest="PysteEnabled", default=True, help="Pyste is only run if this option is supplied to the build command")
	
	gr = opt.parser.add_option_group("pyprop configuration options")
	opt.pyprop_configuration_options = gr
	gr.add_option("--disable-trilinos", action="store_false", dest="PypropUseTrilinos", default=True, help="Disables integration with trilinos completely.")
	gr.add_option("--disable-trilinos-tpetra", action="store_false", dest="PypropUseTrilinosTpetra", default=True, help="Disables trilinos/tpetra integration. Only applies if trilinos is enabled")

	#options from pyprop_waf
	opt.tool_options("pyprop_waf", tooldir="./pyprop/pyprop/build")

	#options from subfolders
	opt.sub_options("einpartikkel")


def configure(conf):
	print('-> Configuring')

	conf.env.CPPFLAGS.append("-I./pyprop/pyprop/core")

	#check dependencies for tool compiler_cxx
	#conf.check_tool("compiler_cxx")
	conf.check_tool("pyprop_waf", tooldir="./pyprop/pyprop/build")

	#configure modules
	conf.sub_config("einpartikkel")


def build(bld):
	bld.add_subdirs("einpartikkel") 


