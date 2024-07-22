set -e

trap cleanup EXIT ERR INT

cleanup(){
if [ "$(ls -A ${TmpD})" ]
  then rm -rf  "${TmpD}"
fi
}

TmpD=$(mktemp -d)
git clone https://github.com/ncbi/sra-tools.git $TmpD/sra-tools
g++ -std=c++17 -O3 -fopenmp -I./src/ -I $TmpD/sra-tools/libs/inc   ./src/aligns_to.cpp ./src/reader.cpp -o ./bin/aligns_to -Wreturn-type -msse4.2 -DBMSSE42OPT -lpthread -ldl -pthread
g++ -std=c++11 -O3 -fopenmp ./src/build_index_of_each_file.cpp -o ./bin/build_index_of_each_file -Wreturn-type
g++ -std=c++11 -O3 -fopenmp ./src/merge_db.cpp -o ./bin/merge_db -Wreturn-type
g++ -std=c++11 -O3 -fopenmp ./src/identify_tax_ids.cpp -o ./bin/identify_tax_ids -Wreturn-type
g++ -std=c++11 -O3 -fopenmp ./src/db_tax_id_to_dbs.cpp -o ./bin/db_tax_id_to_dbs -Wreturn-type
g++ -std=c++11 -O3 -fopenmp ./src/sort_dbs.cpp -o ./bin/sort_dbs -Wreturn-type
echo "SUCCESS"
