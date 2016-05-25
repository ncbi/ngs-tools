cmake_minimum_required (VERSION 2.8)

project ( test-ngs-vdb+c++ )

include_directories ( ${CMAKE_SOURCE_DIR} )

add_executable ( test-ngs-vdb test-ngs-vdb.cpp )
add_dependencies ( test-ngs-vdb ngs-vdb-c++ )

target_link_libraries ( test-ngs-vdb
    ngs-vdb
    ${CMAKE_STATIC_LIBRARY_PREFIX}ktst${CMAKE_STATIC_LIBRARY_SUFFIX}
	${SYS_LIBRARIES}
)

# white box test
add_test ( NgsVdbWhiteBox test-ngs-vdb )

