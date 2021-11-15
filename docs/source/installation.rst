Installation of the python package
==================================

Install the dependencies
------------------------
Required :

- yaml-cpp
- hdf5

Useful to run the example notebooks :

- fklab-python-core
- scipy
- numpy
- matplotlib

Install the package with pip
----------------------------

Soon to be released !


Install from the repository
---------------------------

#. **Clone fklab-compressed-decoder code repositories**

    .. code-block:: console

        cd <PATH>
        git clone https://<USER>@bitbucket.org/kloostermannerflab/fklab-compressed-decoder.git


#. **"Install" repositories**

    .. code-block:: console

        cd fklab-compressed-decoder
        python setup.py build_ext --inplace
        pip install -e . --no-deps


    This will first build all C/C++ extensions in place and add a link to the repository in
    Python's module search path.


Install the fklab-python-core (use in the guide later)
______________________________________________________

Some useful modules have been added in the develop branch of the fklab-python-core repository but not yet released.
To be able to use it, the package needs to be added in develop mode and checkout on the develop branch.

#. **Clone fklab-python-core code repositories**

    .. code-block:: console

        cd <PATH>
        git clone https://<USER>@bitbucket.org/kloostermannerflab/fklab-python-core.git


#. **"Install" repositories**


    .. code-block:: console

        cd fklab-python-core
        git checkout develop
        python setup.py build_ext --inplace
        pip install -e . --no-deps