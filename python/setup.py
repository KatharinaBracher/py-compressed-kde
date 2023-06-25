from setuptools import setup, find_packages
from setuptools.extension import Extension

from distutils.command.build_ext import build_ext
from distutils.sysconfig import get_config_var, customize_compiler

import glob
import os
import sys
import subprocess

root_path = os.path.dirname(__file__)

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
    
src_path =  os.path.abspath(os.path.join(root_path, '../src'))
generated_path = os.path.join(src_path, 'generated')
hgf_path = os.path.join(generated_path, 'HighFive')

schema_path = os.path.join(src_path, 'schema.fbs')


if os.path.exists(generated_path):
    import shutil
    shutil.rmtree(generated_path)
    
os.makedirs(generated_path, exist_ok=True)
os.system(f"flatc --cpp -o {generated_path} {schema_path}")
os.system(f"git clone --depth 1 --quiet --branch v2.3.1 -c advice.detachedHead=false https://github.com/BlueBrain/HighFive.git {hgf_path}")

sources = glob.glob(os.path.abspath(os.path.join(root_path, '../pybind/*.cpp'))) + glob.glob(os.path.join(src_path, '*.cpp'))

include_dirs = [
    src_path, 
    generated_path,
    os.path.join(hgf_path, 'include'),
    get_config_var('INCLUDEDIR'),
    get_pybind_include(), get_pybind_include(user=True), # If compiling against Ubuntu's hdf5, add "/usr/include/hdf5/serial"
]

library_dirs = []

if sys.platform.startswith('win'):
    compile_args = ['/std:c++17', '-DH5_BUILT_AS_DYNAMIC_LIB', '/O2']
    include_dirs.append(os.path.join(get_config_var('prefix'), 'Library', 'include'))
    library_dirs.append(os.path.join(get_config_var('prefix'), 'Library', 'lib'))
elif sys.platform.startswith('linux'):
    compile_args = ['-std=c++17', '-O3']
else:
    compile_args = []


extensions = [
    Extension(
        "compressed_kde.compressed_kde",
        sources = sources,
        libraries = ['yaml-cpp', 'hdf5', 'flatbuffers'], # If compiling against Ubuntu's hdf5, use hdf5_serial
        include_dirs = include_dirs,
        library_dirs = library_dirs,
        language = "c++",
        extra_compile_args = compile_args,
    )
]

import re

VERSIONFILE = os.path.join(root_path, "_version.py")
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))

setup(
    name = "py-compressed-kde",
    version = verstr,
    packages = ['compressed_kde', 'compressed_kde.decode'],
    package_dir={
            "compressed_kde": os.path.join(root_path,"compressed_kde"),
            "compressed_kde.decode": os.path.join(root_path,"compressed_kde/decode"),
            
            },
    install_requires=['h5py', 'pyyaml', 'flatbuffers', "pybind11>=2.10"],
    setup_requires=['pybind11>=2.10'],
    ext_modules = extensions,
    cmdclass = {'build_ext': my_build_ext},
    license_files = ( os.path.join(root_path,'../LICENSE.txt'))
)
    
