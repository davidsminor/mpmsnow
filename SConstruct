import VersionControl
VersionControl.setVersion( "IEBuild", "6.6.4" )

import IEBuild
import os
import glob
import IEEnv

name = "mpmsnow"
version = "1.0.0"

sources = glob.glob( "src/*.cpp" )
includePaths = [ "./include", "/software/tools/include/eigen/3.2.0/cent6.x86_64/include/Eigen" ]

prog = IEBuild.Program( ARGUMENTS, name, sources, includePaths )
prog.setCompiler()

prog.addLib("glut")
prog.addLib("GL")
prog.addLib("GLU")

env = prog.finalize()

installTarget = os.path.join( prog.getInstallPrefix(), "apps", name, version, IEEnv.platform(), "bin" )
env.Install( installTarget, name ) 
env.Alias( "install", installTarget )



