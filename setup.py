long_description = """
Python wrapper for demInterFoamheat.

demInterFoamheat is an OpenFOAM CFD solver which can interact with a
discrete element model.
"""

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import os
import string

assert not os.getenv("FOAM_SRC") is None, "Cannot find OpenFOAM install. Did you source ~/OpenFOAM/OpenFOAM-v3.0+/etc/bashrc ?"

assert os.path.isfile(os.getenv("FOAM_APPBIN")+"/interFoam"), "Cannot find OpenFOAM binaries. Did you build OpenFOAM in the default location? "


def get_version_number(module):
    f = open(os.path.join(module, "__init__.py"))
    return string.strip(f.readline().split("=")[1])[1:-1]

try:
    print "building version", get_version_number("pyDemFoamMFheat")
except:
    print "could not find version number in __init__.py"
    raise


ext = [Extension("_pyDemFoamMFheat",
               sources=["_pyDemFoamMFheat.pyx",
                        "demInterFoamheat.C",
                       "demInterFoamBaseheat.C" ],
              include_dirs = [
                  os.getenv("FOAM_SRC")+"/finiteVolume/lnInclude",
                  os.getenv("FOAM_SRC")+"/meshTools/lnInclude",
                  os.getenv("FOAM_SRC")+"/OpenFOAM/lnInclude",#except
                  os.getenv("FOAM_SRC")+"/OSspecific/POSIX/lnInclude",
                 # os.getenv("FOAM_SOLVERS")+"/multiphase/interFoam",#modified
                  os.getenv("FOAM_SRC")+"/transportModels/twoPhaseMixture/lnInclude",
                  os.getenv("FOAM_SRC")+"/transportModels",
                  os.getenv("FOAM_SRC")+"/transportModels/incompressible/lnInclude",
                  os.getenv("FOAM_SRC")+"/transportModels/interfaceProperties/lnInclude",
                  os.getenv("FOAM_SRC")+"/TurbulenceModels/turbulenceModels/lnInclude",             #turbulence
                  os.getenv("FOAM_SRC")+"/TurbulenceModels/incompressible/lnInclude",               #turbulence
                  os.getenv("FOAM_SRC")+"/transportModels/immiscibleIncompressibleTwoPhaseMixture/lnInclude",
                  os.getenv("FOAM_SRC")+"/sampling/lnInclude",
                  ],
              extra_compile_args= ["-m64", "-Dlinux64", "-DWM_ARCH_OPTION=64",
                                   "-DWM_DP", "-DWM_LABEL_SIZE=32", "-Wall",
                                   "-Wextra", "-Wold-style-cast",
                                   "-Wnon-virtual-dtor",
                                   "-Wno-unused-parameter",
                                   "-Wno-invalid-offsetof", "-O3",
                                   "-DNoRepository", "-ftemplate-depth-100"],
              extra_link_args=["-Xlinker", "--add-needed",
                                  "-Xlinker", "--no-as-needed"],
              libraries = ["immiscibleIncompressibleTwoPhaseMixture","turbulenceModels","incompressibleTurbulenceModels","finiteVolume","fvOptions", "meshTools","sampling","dl","m"],

              library_dirs = [os.getenv("FOAM_LIBBIN")],
              language="c++",             # generate C++ code
    )]


setup(
    name = 'pyDemFoamMFheat',
    packages = ["pyDemFoamMFheat"], # this must be the same as the name above
    version = get_version_number("pyDemFoamMFheat"),
    description = "Python wrapper for Dem OpenFoam solvers.",
    long_description = long_description,
    author = 'Jiayu Lin',
    requires = ['numpy'],
    author_email = 'Jiayu_Lin@tju.edu.cn',
    url = "https://github.com/jkfurtney/PFC3D_OpenFOAM",
    keywords = 'OpenFOAM,CFD,icoFoam,simpleFoam,PFC3D,PFC,DEM'.split(","),
    license          = "BSD",
    classifiers = [
        'Programming Language :: Python :: 2',
        "Development Status :: 5 - Production/Stable",
        "License :: OSI Approved :: BSD License",
        'Topic :: Scientific/Engineering :: Interface Engine/Protocol Translator',
        "Intended Audience :: Science/Research"
    ],
    ext_modules = cythonize(ext, language="c++"))
