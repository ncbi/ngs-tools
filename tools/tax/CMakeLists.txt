if(CMAKE_VERSION VERSION_LESS "3.1")
   if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      set(CMAKE_CXX_FLAGS "-std=gnu++17 ${CMAKE_CXX_FLAGS}")
   endif()
else()
   set(CMAKE_CXX_STANDARD 17)
endif()

FIND_PACKAGE(OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

add_library(ReaderLib STATIC src/reader.cpp src/reader.h)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)


add_executable(aligns_to src/aligns_to.cpp)
target_include_directories(aligns_to PUBLIC src/)
if (UNIX)
target_compile_options(aligns_to PUBLIC -msse4.2 -DBMSSE42OPT)
endif()
target_link_libraries(aligns_to PRIVATE ReaderLib Threads::Threads)

add_executable(tax_collator src/tax_collator.cpp)
target_include_directories(tax_collator PUBLIC src/)
if (UNIX)
target_compile_options(tax_collator PUBLIC -msse4.2 -DBMSSE42OPT)
endif()
target_link_libraries(tax_collator PRIVATE Threads::Threads)


links_and_install_subdir(aligns_to tax)

add_executable(dump_kmers src/dump_kmers.cpp)
target_link_libraries(dump_kmers PRIVATE ReaderLib)
links_and_install_subdir(dump_kmers tax)

add_executable(or_db src/or_db.cpp)
target_link_libraries(or_db PRIVATE ReaderLib)
links_and_install_subdir(or_db tax)

add_executable(and_db src/and_db.cpp)
target_link_libraries(and_db PRIVATE ReaderLib)
links_and_install_subdir(and_db tax)

add_executable(dbs_to_db src/dbs_to_db.cpp)
target_link_libraries(dbs_to_db PRIVATE ReaderLib)
links_and_install_subdir(dbs_to_db tax)

add_executable(sam_filter src/sam_filter.cpp)
target_link_libraries(sam_filter PRIVATE ReaderLib)
links_and_install_subdir(sam_filter tax)

add_executable(subtract_db src/subtract_db.cpp)
target_link_libraries(subtract_db PRIVATE ReaderLib)
links_and_install_subdir(subtract_db tax)

add_executable(subtract_dbs src/subtract_dbs.cpp)
target_link_libraries(subtract_dbs PRIVATE ReaderLib)
links_and_install_subdir(subtract_dbs tax)

add_executable(build_index src/build_index.cpp)
target_link_libraries(build_index PRIVATE ReaderLib)
links_and_install_subdir(build_index tax)

add_executable(build_index_of_each_file src/build_index_of_each_file.cpp)
target_link_libraries(build_index_of_each_file PRIVATE ReaderLib)
links_and_install_subdir(build_index_of_each_file tax)

add_executable(merge_db src/merge_db.cpp)
target_link_libraries(merge_db PRIVATE ReaderLib)
links_and_install_subdir(merge_db tax)

add_executable(merge_tax_ids src/merge_tax_ids.cpp)
target_link_libraries(merge_tax_ids PRIVATE ReaderLib)
links_and_install_subdir(merge_tax_ids tax)

add_executable(merge_kingdoms src/merge_kingdoms.cpp)
target_link_libraries(merge_kingdoms PRIVATE ReaderLib)
links_and_install_subdir(merge_kingdoms tax)

add_executable(build_index_multi src/build_index_multi.cpp)
target_link_libraries(build_index_multi PRIVATE ReaderLib)
links_and_install_subdir(build_index_multi tax)

add_executable(db_fasta_to_bin src/db_fasta_to_bin.cpp)
target_link_libraries(db_fasta_to_bin PRIVATE ReaderLib)
links_and_install_subdir(db_fasta_to_bin tax)

add_executable(db_tax_id_to_dbs src/db_tax_id_to_dbs.cpp)
target_link_libraries(db_tax_id_to_dbs PRIVATE ReaderLib)
links_and_install_subdir(db_tax_id_to_dbs tax)

add_executable(db_to_dbs src/db_to_dbs.cpp)
target_link_libraries(db_to_dbs PRIVATE ReaderLib)
links_and_install_subdir(db_to_dbs tax)

add_executable(identify_tax_ids src/identify_tax_ids.cpp)
target_link_libraries(identify_tax_ids PRIVATE ReaderLib)
links_and_install_subdir(identify_tax_ids tax)

add_executable(db_fasta_to_bin_multi src/db_fasta_to_bin_multi.cpp)
target_link_libraries(db_fasta_to_bin_multi PRIVATE ReaderLib)
links_and_install_subdir(db_fasta_to_bin_multi tax)

add_executable(filter_db src/filter_db.cpp)
target_link_libraries(filter_db PRIVATE ReaderLib)
links_and_install_subdir(filter_db tax)

add_executable(filter_dbs src/filter_dbs.cpp)
target_link_libraries(filter_dbs PRIVATE ReaderLib)
links_and_install_subdir(filter_dbs tax)

add_executable(filter_db_multi src/filter_db_multi.cpp)
target_link_libraries(filter_db_multi PRIVATE ReaderLib)
links_and_install_subdir(filter_db_multi tax)

add_executable(fasta_contamination src/fasta_contamination.cpp)
target_link_libraries(fasta_contamination PRIVATE ReaderLib)
links_and_install_subdir(fasta_contamination tax)

add_executable(fasta_contamination_multi src/fasta_contamination_multi.cpp)
target_link_libraries(fasta_contamination_multi PRIVATE ReaderLib)
links_and_install_subdir(fasta_contamination_multi tax)

add_executable(find_closest_profile_linear src/find_closest_profile_linear.cpp)
target_link_libraries(find_closest_profile_linear PRIVATE ReaderLib)
links_and_install_subdir(find_closest_profile_linear tax)

add_executable(sort_dbs src/sort_dbs.cpp)
target_link_libraries(sort_dbs PRIVATE ReaderLib)
links_and_install_subdir(sort_dbs tax)

add_executable(print_dbs src/print_dbs.cpp)
target_link_libraries(print_dbs PRIVATE ReaderLib)
links_and_install_subdir(print_dbs tax)

include_directories(${CMAKE_SOURCE_DIR})

install(TARGETS aligns_to dump_kmers build_index build_index_of_each_file merge_db merge_tax_ids merge_kingdoms build_index_multi db_to_dbs db_tax_id_to_dbs identify_tax_ids db_fasta_to_bin db_fasta_to_bin_multi filter_db filter_dbs filter_db_multi
                  fasta_contamination fasta_contamination_multi find_closest_profile_linear
                  print_dbs sort_dbs and_db or_db subtract_db subtract_dbs dbs_to_db sam_filter
          RUNTIME DESTINATION bin/tax)

if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    add_executable(get_profile src/get_profile)
    target_link_libraries(get_profile PRIVATE ReaderLib)
    install(TARGETS get_profile RUNTIME DESTINATION bin/tax)
endif()

add_subdirectory(src/tests)
