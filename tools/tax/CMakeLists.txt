if (CMAKE_VERSION VERSION_LESS "3.1")
   if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      set (CMAKE_CXX_FLAGS "-std=gnu++11 ${CMAKE_CXX_FLAGS}")
   endif ()
else ()
   set (CMAKE_CXX_STANDARD 11)
endif ()

FIND_PACKAGE ( OpenMP REQUIRED )
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
set ( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}" )

set ( SHARED_OBJECTS src/reader.cpp )

add_executable ( aligns_to                      src/aligns_to.cpp ${SHARED_OBJECTS})
add_executable ( build_index                    src/build_index.cpp ${SHARED_OBJECTS})
add_executable ( build_index_multi              src/build_index_multi.cpp ${SHARED_OBJECTS})
add_executable ( check_index                    src/check_index.cpp ${SHARED_OBJECTS})
add_executable ( check_index_multi              src/check_index_multi.cpp ${SHARED_OBJECTS})
add_executable ( db_fasta_to_bin                src/db_fasta_to_bin.cpp ${SHARED_OBJECTS})
add_executable ( db_fasta_to_bin_multi          src/db_fasta_to_bin_multi.cpp ${SHARED_OBJECTS})
add_executable ( filter_db                      src/filter_db.cpp ${SHARED_OBJECTS})
add_executable ( filter_db_multi                src/filter_db_multi.cpp ${SHARED_OBJECTS})
add_executable ( contig_builder                 src/contig_builder.cpp ${SHARED_OBJECTS})
add_executable ( fasta_contamination            src/fasta_contamination.cpp ${SHARED_OBJECTS})
add_executable ( find_closest_profile_linear    src/find_closest_profile_linear.cpp ${SHARED_OBJECTS})
add_executable ( contig_connectivity            src/contig_connectivity.cpp ${SHARED_OBJECTS})
add_executable ( sort_dbs                       src/sort_dbs.cpp ${SHARED_OBJECTS})

include_directories ( ${CMAKE_SOURCE_DIR} )

target_link_libraries ( aligns_to ${SYS_LIBRARIES} )
target_link_libraries ( build_index ${SYS_LIBRARIES} )
target_link_libraries ( build_index_multi ${SYS_LIBRARIES} )
target_link_libraries ( check_index ${SYS_LIBRARIES} )
target_link_libraries ( check_index_multi ${SYS_LIBRARIES} )
target_link_libraries ( db_fasta_to_bin ${SYS_LIBRARIES} )
target_link_libraries ( db_fasta_to_bin_multi ${SYS_LIBRARIES} )
target_link_libraries ( filter_db ${SYS_LIBRARIES} )
target_link_libraries ( filter_db_multi ${SYS_LIBRARIES} )
target_link_libraries ( contig_builder ${SYS_LIBRARIES} )
target_link_libraries ( fasta_contamination ${SYS_LIBRARIES} )
target_link_libraries ( find_closest_profile_linear ${SYS_LIBRARIES} )
target_link_libraries ( contig_connectivity ${SYS_LIBRARIES} )
target_link_libraries ( sort_dbs ${SYS_LIBRARIES} )

install ( TARGETS aligns_to build_index build_index_multi check_index check_index_multi db_fasta_to_bin db_fasta_to_bin_multi filter_db filter_db_multi
                  contig_builder fasta_contamination find_closest_profile_linear
                  contig_connectivity sort_dbs
          RUNTIME DESTINATION bin/tax )

if (NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    add_executable ( get_profile src/get_profile ${SHARED_OBJECTS})
    target_link_libraries ( get_profile ${SYS_LIBRARIES} )
    install ( TARGETS get_profile RUNTIME DESTINATION bin/tax )
endif()

add_subdirectory ( src/tests )
