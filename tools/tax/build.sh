sudo apt update
sudo apt install build-essentials
sudo apt install g++
sudo apt install git
sudo apt install make
sudo apt install cmake
sudo apt install default-jdk
sudo apt install ant

git clone https://github.com/ncbi/ngs-tools.git --branch tax
git clone https://github.com/ncbi/ncbi-vdb.git
git clone https://github.com/ncbi/ngs.git
cd ./ngs/
./configure --without-debug
make
cd ..

cd ./ncbi-vdb/
./configure --without-debug
make
cd ..

cd ./ngs-tools/
./configure --without-debug
make
cd ../

cp ~/ncbi-outdir/ngs-tools/linux/gcc/x86_64/rel/bin/* ./ngs-tools/tools/tax/bin
