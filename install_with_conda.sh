conda create -n dp-test python=3.7
conda activate dp-test

export CONDA_OVERRIDE_CUDA=11.3

# install libtensorflow_cc gpu version (Note that only libtensorflow_cc in the deepmodeling channel enables GPU. The libtensorflow_cc library in the conda-forge channel can only be used on GPU.)
conda install -c deepmodeling libtensorflow_cc=*=cuda11*

# install gcc compiler (Note that only pkgs/main works. The conda-forge does not work due to GLIBC version issues. Do not include openmpi-mpicc and openmpi-mpicxx here, or only the cos6 gcc will be installed, and then DeePMD-kitv2.0.3 will not be successfully built. 10/29/2021)
conda install -c pkgs/main gcc_linux-64 gxx_linux-64 openmpi

# install cudatoolkit-dev (Note that cudatoolkit package does not include an nvcc compiler. To have nvcc installed, cudatoolkit-dev is needed.)
conda install -c conda-forge cudatoolkit-dev cudnn
conda install numpy=1.19 cmake fftw
# export LIBRARY_PATH=$CONDA_PREFIX/lib:$LIBRARY_PATH
# export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

##  TF2.6.0 would require C++14 standards, so we need to do some tricks for the compilation of plumed and lammps
pip install tensorflow==2.6.0


# check source/lmp/env.sh.in
# make sure NNP_INC=" -std=c++@CMAKE_CXX_STANDARD@ -D@prec_def@ @TTM_DEF@ -DLAMMPS_VERSION_NUMBER=@LAMMPS_VERSION_NUMBER@ -I$TF_INCLUDE_DIRS -I$DEEPMD_ROOT/include/ "
# instead of NNP_INC=" -std=c++11 -D@prec_def@ @TTM_DEF@ -DLAMMPS_VERSION_NUMBER=@LAMMPS_VERSION_NUMBER@ -I$TF_INCLUDE_DIRS -I$DEEPMD_ROOT/include/ "
# This is why we checkout the devel branch, in the future this will be in master.

git clone -b devel --recursive https://github.com/deepmodeling/deepmd-kit.git deepmd-kit
git clone -b stable --recursive https://github.com/lammps/lammps.git lammps-stable
## lammps are only compatible with plumed v2.4/v2.5/v2.6
git clone -b v2.6.3 --recursive https://github.com/plumed/plumed2.git plumed-v2.6
git clone https://github.com/y1xiaoc/deepmd-plumed.git dp-plmd


tensorflow_root=$CONDA_PREFIX # this cannot be changed
deepmd_root=$CONDA_PREFIX # you can change this
plumed_root=$CONDA_PREFIX # you can change this
lammps_root=$CONDA_PREFIX # you can change this


# dp python
cd deepmd-kit
# do not set the following if you do not use cuda
export DP_VARIANT=cuda
# if there is package conflicts or missing. check duplicates in ~/.local which has higher priority than current env.
pip install  .
# pip install -e .


# dp c++
cd source
mkdir build 
cd build
# remove the last cuda flag if you are not using cuda
cmake ../ -D CMAKE_INSTALL_PREFIX=$CONDA_PREFIX -D TENSORFLOW_ROOT=$CONDA_PREFIX -D USE_CUDA_TOOLKIT=true -D CUDA_TOOLKIT_ROOT_DIR=$CONDA_PREFIX
# cmake ../ -DTENSORFLOW_ROOT=$tensorflow_root -DCMAKE_CXX_FLAGS="-std=c++14" -DCMAKE_INSTALL_PREFIX=$deepmd_root -DUSE_CUDA_TOOLKIT=true -DCUDA_TOOLKIT_ROOT_DIR=$CONDA_PREFIX

make -j 20
make install
make lammps
cd ../../../


# below are that you want to use deepmd model with plumed
cp -r dp-plmd/deepmd/ plumed-v2.6/src/
cd plumed-v2.6
sed -i 's/std=c++11/std=c++14/g' configure
# you may add other modules or flags, and remove "-ldeepmd_op_cuda" if you do not need cuda
./configure --prefix=$plumed_root --enable-modules=dimred:crystallization:deepmd:ves CXX=mpic++ CXXFLAGS="-DHIGH_PREC -O3 -I${CONDA_PREFIX}/include -L${CONDA_PREFIX}/lib -Wl,--no-as-needed -ldeepmd_op -ldeepmd_op_cuda -lrt -ldeepmd -ldeepmd_cc -ltensorflow_cc -ltensorflow_framework -Wl,-rpath=${CONDA_PREFIX}/lib -Wl,-rpath-link=${CONDA_PREFIX}/lib"
# or if you are trying to install things into different places 
# ./configure --prefix=$plumed_root --enable-modules=dimred:crystallization:deepmd CXX=mpic++ CXXFLAGS="-DHIGH_PREC -O3 -I${deepmd_root}/include -I${tensorflow_root}/include -L${deepmd_root}/lib -L${tensorflow_root}/lib -Wl,--no-as-needed -lrt -ldeepmd_op -ldeepmd_op_cuda -ldeepmd -ldeepmd_cc -ltensorflow_cc -ltensorflow_framework -Wl,-rpath=${tensorflow_root}/lib -Wl,-rpath-link=${tensorflow_root}/lib -Wl,-rpath=${deepmd_root}/lib -Wl,-rpath-link=${deepmd_root}/lib"
make -j 20
make install
cd ../


# lammps 
cd lammps-stable
cp -r ../deepmd-kit/source/build/USER-DEEPMD/ src/ 
sed -i 's/set(STANDARD_PACKAGES/set(STANDARD_PACKAGES USER-DEEPMD/g' cmake/CMakeLists.txt
sed -i 's/set(CMAKE_CXX_STANDARD 11)/set(CMAKE_CXX_STANDARD 14)/g' cmake/CMakeLists.txt

mkdir build
cd build
# add the package you want to use here, remove plumed related if you do not use plumed
ARGS="-D PKG_MANYBODY=yes -D PKG_DIPOLE=ON -D PKG_KSPACE=ON -D PKG_MISC=ON -D PKG_MOLECULE=ON -D PKG_RIGID=ON -D PKG_OPT=ON -D PKG_PLUMED=yes -D DOWNLOAD_PLUMED=no -D PLUMED_MODE=shared"
# remove "-ldeepmd_op_cuda" if you do not need cuda
cmake -D CMAKE_BUILD_TYPE=Release $ARGS -D PKG_USER-DEEPMD=yes -D FFT=FFTW3 -D CMAKE_INSTALL_PREFIX=${lammps_root} -D CMAKE_CXX_FLAGS="-DHIGH_PREC -I${CONDA_PREFIX}/include -L${CONDA_PREFIX}/lib -Wl,--no-as-needed -lrt -ldeepmd_op -ldeepmd_op_cuda -ldeepmd -ldeepmd_cc -ltensorflow_cc -ltensorflow_framework -Wl,-rpath=${CONDA_PREFIX}/lib -Wl,-rpath-link=${CONDA_PREFIX}/lib -DLAMMPS_VERSION_NUMBER=20210920" ../cmake
make -j 20
make install
cd ../../


# you may need to set the mpi cuda support to true before you run mpi+cuda programs
