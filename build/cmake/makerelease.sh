mkdir -p Release
cd Release
cmake -DCMAKE_INSTALL_PREFIX:PATH=~/install/ ../../..
make
