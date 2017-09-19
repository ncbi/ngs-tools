include_directories ( ${CMAKE_SOURCE_DIR} )

add_executable ( test-ngs-vdb test-ngs-vdb.cpp )
add_dependencies ( test-ngs-vdb ngs-vdb )

target_link_libraries ( test-ngs-vdb
    ngs-vdb
    tk-version
    ${CMAKE_STATIC_LIBRARY_PREFIX}ktst${CMAKE_STATIC_LIBRARY_SUFFIX}
	${SYS_LIBRARIES}
)

# white box test
add_test ( NAME NgsVdbWhiteBox COMMAND test-ngs-vdb )

