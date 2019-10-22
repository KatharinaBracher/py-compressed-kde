from setuptools import setup
from setuptools.extension import Extension

from distutils.command.build_ext import build_ext
from distutils.sysconfig import get_config_var, customize_compiler

import glob
import os


# use custom build_ext class that removes the -Wstrict-prototypes compiler flag
# this flag is not supported for C++ and results in warnings
class my_build_ext(build_ext):
    def build_extensions(self):
        customize_compiler(self.compiler)
        try:
            self.compiler.compiler_so.remove("-Wstrict-prototypes")
        except (AttributeError, ValueError):
            pass
        build_ext.build_extensions(self)


# helper class to get include directories for pybind11
class get_pybind_include(object):
    def __init__(self, user=False):
        self.user = user
    
    def __iter__(self):
        for k in str(self):
            yield k
    
    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)
  
sources = glob.glob(os.path.abspath('../pybind/*.cpp')) + glob.glob(os.path.abspath('../src/*.cpp'))

extensions = [
    Extension(
        "fklab.decode.compressed_kde",
        sources = sources,
        libraries = ['yaml-cpp', 'hdf5'],
        include_dirs = [os.path.abspath('../src'), os.path.abspath('../ext/HighFive-1.4/include'), get_config_var('INCLUDEDIR'), get_pybind_include(), get_pybind_include(user=True)],
        language = "c++",
        extra_compile_args = ['-std=c++14', '-O3'],
    )
]

import re

VERSIONFILE = "fklab/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))


setup(
    name = "py-compressed-decoder",
    version = verstr,
    packages = ['fklab.decode'],
    install_requires=['hdf5', 'yaml-cpp', 'boost'],
    ext_modules = extensions,
    cmdclass = {'build_ext': my_build_ext},
)
    
