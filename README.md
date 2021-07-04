# DeepMD Plumed Module

This module provides interface to use DeePMD models (deep potential and deep tensor) as plumed CVs. 

## Installation

This module is to be patched to plumed source code and compiled together. 
To use it, first make sure deepmd and tensorflow has been installed correctly. Then copy the `deepmd` folder into the `src` folder of plumed directory. 
Then see the following script:
```bash
# do not set the following if you do not use cuda
DEEPMD_CUDA_LINK="-ldeepmd_op_cuda"

# at plumed source folder
./configure --prefix=$plumed_root --enable-modules=+deepmd CXX=mpic++ CXXFLAGS="-DHIGH_PREC -O3 -I${tensorflow_root}/include -I${deepmd_root}/include -L${tensorflow_root}/lib -L${deepmd_root}/lib -Wl,--no-as-needed -ldeepmd_op ${DEEPMD_CUDA_LINK} -lrt -ldeepmd -ldeepmd_cc -ltensorflow_cc -ltensorflow_framework -Wl,-rpath=${tensorflow_root}/lib -Wl,-rpath=${deepmd_root}/lib -Wl,-rpath-link=${tensorflow_root}/lib -Wl,-rpath-link=${deepmd_root}/lib"
```
Note here we follow the same convension as installing DeePMD, let `plumed_root` be the place you want to install plumed to. 
And `tensorflow_root` and `deepmd_root` should be set properly to their installed path.
After that you can `make` and `make install` as usual.

An example install script that install deepmd, lammps and plumed all into a conda environment can be found [here](./install_with_conda.sh). Note it will use gpu and install cuda related things. You may need to modify it to fit your needs.

## Usage

Three different CVs are provided by this module, all using deepmd models.
- `DEEPPOTENTIAL`: Calculate the potential energy (or any scalar model).
- `DEEPDIPOLE`: Calculate the dipole moment (or any vector).
- `DEEPPOLAR`: Calculate the polarizability (or any rank 2 tensor).

All calculations are giving global outputs, no atomic value is available. 
For dipole there are 3 components `x`, `y` and `z`.
For polar there are 9, `xx`, `xy`, `xz`, `yx`, `yy`, `yz`, `zx`, `zy`, `zz`.

A binary graph file of the model (`MODEL` keyword) and 
a text file (`ATYPE` keyword) specify the type of each atom are needed.
The type file should be a list of integers with the i-th element correspond to the type 
in the deep tensor model of the i-th atom in the system (or specified by you).
The output is scaled by a factor given by `UNIT_CVT`, 
which defaults to 1 for tensor models and convert to plumed energy unit for potential model.

By default will use periodic boundary conditions, which will be handled by 
the deep tensor model automatically. In case `NOPBC` flag is specified, the box
will be ignored and there will be no the pbc handling.

Here's a simple example showing how to use these CV (also the default values of the keywords):
```
dp: DEEPPOTENTIAL MODEL=graph.pb ATYPE=type.raw # UNIT_CVT=96.487 # by default convert to plumed unit
dipole: DEEPDIPOLE MODEL=dipole.pb ATYPE=type.raw UNIT_CVT=1.0
polar: DEEPPOLAR MODEL=polar.pb ATYPE=type.raw UNIT_CVT=1.0
```
