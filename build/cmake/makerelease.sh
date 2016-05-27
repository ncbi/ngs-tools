mkdir -p Release
cd Release
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX:PATH=~/install/ ../../..
make
