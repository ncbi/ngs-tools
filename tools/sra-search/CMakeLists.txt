add_executable ( sra-search main.cpp vdb-search.cpp searchblock.cpp )

include_directories ( ${CMAKE_SOURCE_DIR} )

add_dependencies ( sra-search ngs-vdb )

target_link_libraries ( sra-search
    ngs-vdb
	${SYS_LIBRARIES}
)

install ( TARGETS sra-search RUNTIME DESTINATION bin )

add_subdirectory ( test )
