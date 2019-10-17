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
        libraries = ['yaml-cpp==0.6.2', 'hdf5'],
        include_dirs = [os.path.abspath('../src'), os.path.abspath('../ext/HighFive-1.4/include'), get_config_var('INCLUDEDIR'), get_pybind_include(), get_pybind_include(user=True)],
        language = "c++",
        extra_compile_args = ['-std=c++14', '-O3'],
    )
]

setup(
    name = "py-compressed-decoder",
    version = "0.2.0",
    packages = ['fklab.decode'],
    install_requires=['hdf5', 'yaml-cpp==0.6.2', 'boost'],
    ext_modules = extensions,
    cmdclass = {'build_ext': my_build_ext},
)
    
