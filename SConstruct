import VersionControl
VersionControl.setVersion( "IEBuild", "6.6.4" )

import IEBuild
import os
import glob
import IEEnv

name = "mpmsnow"
version = "1.0.0"

sources = [ "src/main.cpp" ] + glob.glob( "src/MpmSim/*.cpp" )
includePaths = [ "./include", "/software/tools/include/eigen/3.2.0/cent6.x86_64/include/Eigen" ]

prog = IEBuild.Program( ARGUMENTS, name, sources, includePaths )
prog.setCompiler()

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

env = prog.finalize()

installTarget = os.path.join( prog.getInstallPrefix(), "apps", name, version, IEEnv.platform(), "bin" )
env.Install( installTarget, name ) 
env.Alias( "install", installTarget )



