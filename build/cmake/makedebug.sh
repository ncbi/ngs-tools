mkdir -p Debug
cd Debug
cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX:PATH=~/usr/ ../../..
make
