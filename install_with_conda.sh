conda create -n dpdev
source activate dpdev
conda install -c deepmodeling -c conda-forge -c nvidia --override-channels libtensorflow_cc=*=gpu_cuda11* cudatoolkit=11  gxx_linux-64 gcc_linux-64 gfortran_linux-64 openmpi=4 openmpi-mpicc openmpi-mpicxx openmpi-mpifort fftw cmake gtest gmock
# make sure numpy is 1.19 to avoid conflict with tf
conda install -c conda-forge numpy=1.19 scipy ipython jupyter matplotlib
# or if you want mkl
# conda install -c conda-forge numpy=1.19 scipy ipython jupyter matplotlib "libblas=*=*mkl"

pip install tensorflow

git clone -b devel --recursive https://github.com/deepmodeling/deepmd-kit.git deepmd-kit
git clone -b stable --recursive https://github.com/lammps/lammps.git lammps-stable
git clone -b v2.6.3 --recursive https://github.com/plumed/plumed2.git plumed-v2.6
git clone https://github.com/y1xiaoc/deepmd-plumed.git dp-plmd

# restart shell here to avoid strange errors

# now we setup all paths needed
# you may want to also set those in your conda env's activation.d
tensorflow_root=$CONDA_PREFIX # this cannot be changed
deepmd_root=$CONDA_PREFIX # you can change this
plumed_root=$CONDA_PREFIX # you can change this
lammps_root=$CONDA_PREFIX # you can change this


# dp python
cd deepmd-kit
# do not set the following if you do not use cuda
export DP_VARIANT=cuda
# if you just want a simple install
pip install .
# or add -e to install in dev mode (file change take effect immediately)
# pip install -e .[test]


# dp c++
cd source
mkdir build 
cd build
# remove the last cuda flag if you are not using cuda
cmake -DTENSORFLOW_ROOT=$tensorflow_root -DCMAKE_INSTALL_PREFIX=$deepmd_root -DUSE_CUDA_TOOLKIT=true ..
make -j 10
make install
make lammps
cd ../../../


# below are that you want to use deepmd model with plumed
cp -r dp-plmd/deepmd/ plumed-v2.6/src/
cd plumed-v2.6
# you may add other modules or flags, and remove "-ldeepmd_op_cuda" if you do not need cuda
./configure --prefix=$plumed_root --enable-modules=dimred:crystallization:deepmd CXX=mpic++ CXXFLAGS="-DHIGH_PREC -O3 -I${CONDA_PREFIX}/include -L${CONDA_PREFIX}/lib -Wl,--no-as-needed -ldeepmd_op -ldeepmd_op_cuda -lrt -ldeepmd -ldeepmd_cc -ltensorflow_cc -ltensorflow_framework -Wl,-rpath=${CONDA_PREFIX}/lib -Wl,-rpath-link=${CONDA_PREFIX}/lib"
# or if you are trying to install things into different places 
# ./configure --prefix=$plumed_root --enable-modules=dimred:crystallization:deepmd CXX=mpic++ CXXFLAGS="-DHIGH_PREC -O3 -I${deepmd_root}/include -I${tensorflow_root}/include -L${deepmd_root}/lib -L${tensorflow_root}/lib -Wl,--no-as-needed -lrt -ldeepmd_op -ldeepmd_op_cuda -ldeepmd -ldeepmd_cc -ltensorflow_cc -ltensorflow_framework -Wl,-rpath=${tensorflow_root}/lib -Wl,-rpath-link=${tensorflow_root}/lib -Wl,-rpath=${deepmd_root}/lib -Wl,-rpath-link=${deepmd_root}/lib"
make -j 10
make install
cd ../


# lammps 
cd lammps-stable
cp -r ../deepmd-kit/source/build/USER-DEEPMD/ src/ 
sed -i 's/set(STANDARD_PACKAGES/set(STANDARD_PACKAGES USER-DEEPMD/g' cmake/CMakeLists.txt
mkdir build
cd build
# add the package you want to use here, remove plumed related if you do not use plumed
ARGS="-D PKG_MANYBODY=yes -D PKG_DIPOLE=ON -D PKG_KSPACE=ON -D PKG_MISC=ON -D PKG_MOLECULE=ON -D PKG_RIGID=ON -D PKG_OPT=ON -D PKG_USER-MISC=ON -D PKG_USER-PLUMED=yes -D DOWNLOAD_PLUMED=no -D PLUMED_MODE=shared"
# remove "-ldeepmd_op_cuda" if you do not need cuda
cmake -D CMAKE_BUILD_TYPE=Release $ARGS -D PKG_USER-DEEPMD=yes -D FFT=FFTW3 -D CMAKE_INSTALL_PREFIX=${lammps_root} -D CMAKE_CXX_FLAGS="-DHIGH_PREC -I${CONDA_PREFIX}/include -L${CONDA_PREFIX}/lib -Wl,--no-as-needed -lrt -ldeepmd_op -ldeepmd_op_cuda -ldeepmd -ldeepmd_cc -ltensorflow_cc -ltensorflow_framework -Wl,-rpath=${CONDA_PREFIX}/lib -Wl,-rpath-link=${CONDA_PREFIX}/lib" ../cmake
# or if you are trying to install things into different places
#cmake -D CMAKE_BUILD_TYPE=Release -D BUILD_LIB=on -D BUILD_SHARED_LIBS=on -DCMAKE_INSTALL_LIBDIR=lib $ARGS -D PKG_USER-DEEPMD=ON -D FFT=FFTW3 -D CMAKE_INSTALL_PREFIX=${lammps_root} -D CMAKE_CXX_FLAGS="-DHIGH_PREC -I${deepmd_root}/include -I${tensorflow_root}/include -I${plumed_root}/includ -L${deepmd_root}/lib -L${tensorflow_root}/lib -L${plumed_root}/lib -Wl,--no-as-needed -lrt -ldeepmd_op -ldeepmd_op_cuda -ldeepmd -ldeepmd_cc -ltensorflow_cc -ltensorflow_framework -Wl,-rpath=${tensorflow_root}/lib -Wl,-rpath-link=${tensorflow_root}/lib -Wl,-rpath=${deepmd_root}/lib -Wl,-rpath=${plumed_root}/lib -Wl,-rpath=${tensorflow_root}/lib" ../cmake
make -j 10
make install
cd ../../


# you may need to set the mpi cuda support to true before you run mpi+cuda programs
# export OMPI_MCA_opal_cuda_support=true
