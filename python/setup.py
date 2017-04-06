class get_pybind_include(object):
    def __init__(self, user=False):
        self.user = user
    
    def __iter__(self):
        for k in str(self):
            yield k
    
    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)
    
    


def configuration(parent_package='', top_path=None):
    
    from numpy.distutils.misc_util import Configuration
    from distutils.sysconfig import get_config_var
    
    config = Configuration('fklab/decode', parent_package, top_path)
    
    config.add_extension( 'compressed_kde',
                          sources = ['../pybind/*.cpp', '../src/*.cpp'],
                          libraries = ['yaml-cpp'],
                          include_dirs = ['../src', get_config_var('INCLUDEDIR'), get_pybind_include(), get_pybind_include(user=True)],
                          language = "c++",
                          extra_compile_args = ['-std=c++11', '-O3'])
    
    return config
    

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
