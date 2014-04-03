import VersionControl
VersionControl.setVersion( "IEBuild", "6.6.4" )

import IEBuild
import os
import glob
import IEEnv

name = "mpmsnow"
version = "1.0.0"


# build shared sim library:
lib = IEBuild.SharedLibrary(
	ARGUMENTS,
	'MpmSim',
	version = None,
	source = glob.glob( "src/MpmSim/*.cpp" ),
	includePath = [
		"include",
		"/software/tools/include/eigen/3.2.0/cent6.x86_64/include/Eigen"
	]
)

lib.setCompiler()
lib.setCortexMajorVersion()
lib.addLib("GL")
lib.addLib( "tbb" )

lib.finalize()



# build test program:
prog = IEBuild.Program(
	ARGUMENTS,
	"mpmsnow",
	[ "src/main.cpp" ],
	[
		"./include",
		"/software/tools/include/eigen/3.2.0/cent6.x86_64/include/Eigen"
	],
)

prog.setCompiler()

prog.addLib("MpmSim")
prog.addLib("glut")
prog.addLib("GL")
prog.addLib("GLU")

prog.setCortexMajorVersion()
prog.addCortexLib( 'IECore' )
prog.addBoostLib( "regex" )
prog.addBoostLib( "system" )
prog.addBoostLib( "filesystem" )
prog.addBoostLib( "iostreams" )
prog.addBoostLib( "thread" )
prog.addBoostLib( "signals" )
prog.addBoostLib( "wave" )
prog.addLib( "tbb" )

prog.addOpenEXRLib( "Iex" )
prog.addOpenEXRLib( "IlmImf" )
prog.addOpenEXRLib( "IlmThread" )
prog.addOpenEXRLib( "Half" )
prog.addOpenEXRLib( "Imath" )
prog.addDefine( 'HAVE_CORTEX', "1" )

prog.finalize()



# build houdini plugin:
plugin = IEBuild.HoudiniPlugin(
	ARGUMENTS,
	'ieMPMSim',
	version,
	source = glob.glob( "src/houdiniPlugin/*.cpp" ),
	includePath = 
	[
		"./include",
		"/software/tools/include/eigen/3.2.0/cent6.x86_64/include/Eigen"
	]
)

## \todo: should these two statements be default for HoudiniPlugin?
plugin.addHoudiniToolIncludes()
plugin.setWarningsAsErrors( False )
plugin.addLibPaths(["."])
plugin.addLib("MpmSim")
plugin.addLib( "openvdb_sesi" )
plugin.addDefine( "OPENVDB_ENABLED" )
plugin.finalize()



