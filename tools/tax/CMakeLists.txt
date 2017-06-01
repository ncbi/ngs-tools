set ( CMAKE_CXX_STANDARD 11 )
FIND_PACKAGE ( OpenMP REQUIRED )
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
set ( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}" )

set ( SHARED_OBJECTS src/reader.cpp )

add_executable ( aligns_to src/aligns_to ${SHARED_OBJECTS})
add_executable ( build_index src/build_index ${SHARED_OBJECTS})
add_executable ( check_index src/check_index ${SHARED_OBJECTS})
add_executable ( db_fasta_to_bin src/db_fasta_to_bin ${SHARED_OBJECTS})
add_executable ( filter_db src/filter_db ${SHARED_OBJECTS})
add_executable ( sort_dbs src/sort_dbs ${SHARED_OBJECTS})

include_directories ( ${CMAKE_SOURCE_DIR} )

target_link_libraries ( aligns_to ${SYS_LIBRARIES} )
target_link_libraries ( build_index ${SYS_LIBRARIES} )
target_link_libraries ( check_index ${SYS_LIBRARIES} )
target_link_libraries ( db_fasta_to_bin ${SYS_LIBRARIES} )
target_link_libraries ( filter_db ${SYS_LIBRARIES} )
target_link_libraries ( sort_dbs ${SYS_LIBRARIES} )

install ( TARGETS aligns_to build_index check_index db_fasta_to_bin filter_db sort_dbs RUNTIME DESTINATION bin/tax )

add_subdirectory ( src/tests )
