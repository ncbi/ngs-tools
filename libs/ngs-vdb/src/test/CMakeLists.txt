include_directories ( ${CMAKE_SOURCE_DIR} )
message(CMAKE_CXX_FLAGS(test-ngs-vdb) = ${CMAKE_CXX_FLAGS})
message(CMAKE_CXX_STANDARD(test-ngs-vdb) = ${CMAKE_CXX_STANDARD})
add_executable ( test-ngs-vdb test-ngs-vdb.cpp )
add_dependencies ( test-ngs-vdb ngs-vdb )

target_link_libraries ( test-ngs-vdb
    ngs-vdb
    tk-version
    ${CMAKE_STATIC_LIBRARY_PREFIX}ktst${CMAKE_STATIC_LIBRARY_SUFFIX}
	${SYS_LIBRARIES}
)

# white box test
add_test ( NAME SlowTest_NgsVdbWhiteBox COMMAND test-ngs-vdb )

