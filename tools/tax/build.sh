sudo apt install git
sudo apt install make 
sudo apt install cmake
#sudo apt install default-jdk
git clone https://github.com/ncbi/ngs-tools.git --branch tax
git clone https://github.com/ncbi/ncbi-vdb.git
git clone https://github.com/ncbi/ngs.git
./ngs/configure --build-prefix=~/ncbi-outdir --without-debug
cd ./ngs/
make
cd ..

./ncbi-vdb/configure --build-prefix=~/ncbi-outdir --without-debug --with-ngs-sdk-prefix=~/ncbi-outdir
cd ./ncbi-vdb/
make
cd ..

./ngs-tools/configure '--without-debug' '--with-ngs-sdk-prefix=~/ncbi-outdir' '--with-ncbi-vdb-build=~/ncbi-outdir'
cd ./ngs-tools/tools/tax
make 
cd ..

cp ~/ncbi-outdir/ngs-tools/linux/gcc/x86_64/rel/bin/* ./ngs-tools/tools/tax/bin
