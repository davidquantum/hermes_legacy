Linux
=====

Download and compilation
~~~~~~~~~~~~~~~~~~~~~~~~

If you are using a Debian-based system, install the required libraries first:: 

    apt-get install cmake g++ gfortran freeglut3-dev libsuitesparse-dev libglew1.5-dev python-dev python-numpy python-scipy cython python-matplotlib

(Note: cmake has to be at least version 2.6 or later, matplotlib has to be at
least 0.98.5.2 or higher.)

There are multiple ways to obtain the source code. The simplest one 
is to clone the Git repository either from the hpfem.org server::

    git clone http://git.hpfem.org/git/hermes.git

or from Github::
  
    git clone git://github.com/hpfem/hermes.git

These two repositories are synchronized. For more advanced users we recommend to 
create a free account at `Github <http://github.com>`_ (if you do not have one yet),
fork the `Hermes repository <http://github.com/hpfem/hermes>`_, and then clone your 
Github copy of Hermes to your local computer. This will establish links between
your local copy and the master repository, and you'll become part of the Hermes 
network at Github.

Once you have a local copy of the Hermes repository on your computer, change dir 
to hermes/. There you will find a CMakeLists.txt file that contains the lines
::

    # Testing.
    set(WITH_TESTS                YES)
    
    # Optional parts of the library:
    set(WITH_H1D                  YES)
        # Release and debug versions.
        set(H1D_DEBUG             YES)
        ...
  
    set(WITH_H2D                  YES)
        # Real/complex version of the library (you can build either one or both):
        set(H2D_REAL              YES)
        set(H2D_COMPLEX           YES)
        # Release and debug versions.
        set(H2D_DEBUG             YES)
        ...
  
    set(WITH_H3D                  YES)
        # Real/complex version of the library.
        set(H3D_REAL              YES)
        ...
    set(WITH_SUPERLU            NO)
    ...


Create a file called "CMake.vars" where you set all 
these variables according to your needs. Examples of CMake.vars files can
be found in the CMakeVars folder.
After that, type::

    cmake .
    make

If you have more than one CPU, you can use "make -jN" where N is
the number of CPUs of your computer.

Tests
~~~~~

To execute all tests, do::

    ctest -jN

where N is the number of your cores. Note that some tests (especially for adaptivity 
algorithms) take a longer time to finish. To just execute the short running tests, type::

    make test-quick

More options
~~~~~~~~~~~~

You can turn on and off various components to build, just create the CMake.vars
file and add the following::

    set(WITH_EXAMPLES NO)
    set(WITH_PYTHON YES)

(and any other option that you would like to change, see CMakeLists.txt for the
whole list).

You can also easily generate it from a script (e.g. a debian/rules file) by::

    python -c 'print "set(H2D_COMPLEX no)\nset(WITH_EXAMPLES no)\nset(WITH_TUTORIAL no)\nset(WITH_PYTHON yes)\nset(WITH_GLUT no)\nset(WITH_UTIL no)"' > CMake.vars

If you are on OS X, you have to disable GLUT as the glut library is not easily
installable on OS X. To do so, just put the following line into your
CMake.vars::

    set(WITH_GLUT NO)

Debugging with Eclipse
~~~~~~~~~~~~~~~~~~~~~~

To use eclipse as debugger, in the root folder of the project::

    mkdir eclipse_build
    cd eclipse_build
    cmake -G"Eclipse CDT4 - Unix Makefiles" -D CMAKE_BUILD_TYPE=Debug ../

In Eclipse:

    - Import project using Menu File->Import
    - Select General->Existing projects into workspace:
    - Browse where your build tree is and select the root build tree directory. 
    - Keep "Copy projects into workspace" unchecked.


Install Hermes
~~~~~~~~~~~~~~

::

    cmake -DCMAKE_INSTALL_PREFIX=~/usr .
    make
    make install
