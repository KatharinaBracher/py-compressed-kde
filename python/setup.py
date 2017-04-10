from setuptools import setup
from setuptools.extension import Extension
from distutils.sysconfig import get_config_var
import glob
import os

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
        libraries = ['yaml-cpp'],
        include_dirs = [os.path.abspath('../src'), get_config_var('INCLUDEDIR'), get_pybind_include(), get_pybind_include(user=True)],
        language = "c++",
        extra_compile_args = ['-std=c++11', '-O3'],
    )
]

setup(
    name = "py-compressed-decoder",
    packages = ['fklab.decode'],
    ext_modules = extensions
)
    
