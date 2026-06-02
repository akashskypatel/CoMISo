![CoMISo](Logo.png)

# CoMISo -- Constrained Mixed-Integer Solver
A handy solver for optimizing discrete quadratic energies subject to linear and integer constraints, performing proper elimination of the constraints, while relieving the user of cumbersome re-indexing. The solver has been successfully deployed in high-end geometry processing tasks such as the *Mixed-Integer Quadrangulation* project.

## Requirements
Here is an example of what packages were needed to compile CoMISo on a freshly installed Ubuntu 9.04 system

    sudo apt-get install g++
    sudo apt-get install cmake
    sudo apt-get install libgmm-dev
    sudo apt-get install libboost-dev
    sudo apt-get install libblas-dev
    sudo apt-get install libsuitesparse-dev 
(some other needed libraries such as lapack, are installed as dependencies of the above)

For Windows and Macintosh systems the corresponding packages need to be downloaded and installed.

The cmake build system should enable building the CoMISo library under Windows and Macintosh systems, please let me know if this is (not) the case!

## OpenFlipper requirements
To build OpenFlipper you additionally need to install all the Qt4 packages libqt4-{dev-dbg, dev, network, gui, opengl, opengl-dev, script, scripttools, ...} and also

    sudo apt-get install libglew1.5-dev
    sudo apt-get install glutg3-dev

## Building (standalone)
Assuming CoMISo was unpacked to the directory `SOME_DIRECTORY/CoMISo` (where `SOME_DIRECTORY` should be `/PATH_TO_OPENFLIPPER/libs/CoMISo` for integration with the OpenFlipper framework) the package is built by creating a build directory, using cmake to create the Makefiles and using make to actually build:

    cd /SOME_DIRECTORY/CoMISo/
    mkdir build
    cd build
    cmake ..

(assuming all needed packages are installed and cmake threw no errors...)

    make

The binaries (examples) and the shared library are found under `/SOME_DIRECTORY/CoMISO/build/Build/bin/` and `/SOME_DIRECTORY/CoMISO/build/Build/lib/CoMISo/`.

## Building via setup.py
This repository also includes a `setup.py` wrapper for the native CMake build. It is intended for environments that want to drive the build through Python while still using CoMISo's existing CMake targets.

Using the requested virtual environment on Windows:

```powershell
python setup.py build_native
python setup.py install_native
# build for linux via wsl
python setup.py build_native_linux
python setup.py build_native_linux --wsl-distro Ubuntu
```

By default, `install_native` uses the active Python environment's `sys.prefix` as `CMAKE_INSTALL_PREFIX`, so running it with the above interpreter installs CoMISo into that virtual environment. Extra CMake configure flags can be passed with `--cmake-args="..."` or the `COMISO_CMAKE_ARGS` environment variable.

`build_native_linux` relays the build through WSL. Pass `--wsl-distro` to target a specific distro, or set `COMISO_WSL_DISTRO`. If neither is provided, the first distro returned by `wsl.exe -l -q` is used automatically. Override the Linux-side Python executable with `COMISO_WSL_PYTHON`.

On a native Linux environment, use the direct native commands instead of the WSL relay:

```bash
python3 setup.py build_native
python3 setup.py install_native
```

The native Linux path runs fully inside the local Linux environment and uses a separate default build directory from Windows (`build/python-linux` versus `build/python-windows`).

On a native macOS environment, use the same direct native commands:

```bash
python3 setup.py build_native
python3 setup.py install_native
```

If `cmake` is not already available on `PATH`, install the Python `cmake` package for the active interpreter first, for example:

```bash
python3 -m pip install --user cmake
```

The native macOS path uses its own default build directory (`build/python-macos`).


## Building (for use with OpenFlipper)
Simply extract / checkout the CoMISo directory to the `/PATH_TO_OPENFLIPPER/libs/` directory. The library will be automatically built and you will find the shared library `libCoMISo.so` under the OpenFlipper build directory.
To use the solver in your plugin, add CoMISo to the `CMakeLists.txt` of the plugin and you are set, see *Plugin-HarmonicExample* for an example.

## Using
To use the solver library in your applications have a look at the `/SOME_DIRECTORY/CoMISo/Examples/` and the sample OpenFlipper plugin (*Plugin-HarmonicExample*) downloadable from the CoMISo project homepage.

## Feedback
We appreciate your feedback! Bugs, comments, questions or patches send them to <zimmer@informatik.rwth-aachen.de> or <bommes@informatik.rwth-aachen.de>!
