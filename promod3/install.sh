# #!/bin/bash

# # Update and install prerequisites
# sudo apt update
# sudo apt install -y build-essential libssl-dev g++ python3.9 python3.9-dev python3.9-venv \
#   autotools-dev libicu-dev libbz2-dev zlib1g-dev sqlite3 libsqlite3-dev \
#   libfftw3-dev libtiff5 libtiff-dev libpng-dev

# # Install CMake 3.23.1
# echo "Installing CMake 3.23.1..."
# wget https://github.com/Kitware/CMake/releases/download/v3.23.1/cmake-3.23.1.tar.gz
# tar -zxvf cmake-3.23.1.tar.gz
# cd cmake-3.23.1
# ./bootstrap
# make
# sudo make install
# cd ..
# rm -rf cmake-3.23.1 cmake-3.23.1.tar.gz

# # Install Boost 1.76.0
# echo "Installing Boost 1.76.0..."
# wget https://boostorg.jfrog.io/artifactory/main/release/1.76.0/source/boost_1_76_0.tar.gz
# tar -zxvf boost_1_76_0.tar.gz
# cd boost_1_76_0
# ./bootstrap.sh
# sudo ./b2 install
# cd ..
# rm -rf boost_1_76_0 boost_1_76_0.tar.gz

# # Install Eigen3 3.4.0
# echo "Installing Eigen3 3.4.0..."
# wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
# tar -xvzf eigen-3.4.0.tar.gz
# cd eigen-3.4.0
# mkdir build
# cd build
# cmake ..
# sudo make install
# cd ../..
# rm -rf eigen-3.4.0 eigen-3.4.0.tar.gz



# # Confirmation message
# echo "All packages have been successfully installed."


sudo apt-get install cmake g++ libtiff-dev libeigen3-dev \
             libpng-dev libboost-all-dev \
             libpng-dev libsqlite3-dev

# Install FFTW3 3.3.9 (with single precision)
echo "Installing FFTW3 3.3.9..."
wget http://www.fftw.org/fftw-3.3.9.tar.gz
tar -xvzf fftw-3.3.9.tar.gz
cd fftw-3.3.9
./configure --enable-single CFLAGS="-fPIC"
make
sudo make install
cd ..
rm -rf fftw-3.3.9 fftw-3.3.9.tar.gz

OPEN_MM_ROOT="/opt/conda/lib/python3.10"
# COMPOUND_LIB="/home/jerry/mara/protein_structure_processing/compound_lib/compounds.chemlib"

# Install openstructure
git clone https://git.scicore.unibas.ch/schwede/openstructure.git
cd openstructure
cmake . -DOPTIMIZE=ON -DENABLE_INFO=OFF -DENABLE_MM=ON \
        -DOPEN_MM_LIBRARY=$OPEN_MM_ROOT/lib/libOpenMM.so \
        -DOPEN_MM_INCLUDE_DIR=$OPEN_MM_ROOT/include/ \
        -DOPEN_MM_PLUGIN_DIR=$OPEN_MM_ROOT/lib/plugins \
        # -DCOMPOUND_LIB=$COMPOUND_LIB
make -j4