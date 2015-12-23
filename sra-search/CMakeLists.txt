cmake_minimum_required (VERSION 2.8)

project ( sra-search )

add_executable ( sra-search main.cpp vdb-search.cpp )

add_dependencies ( sra-search ngs ncbi-vdb ) 

target_link_libraries ( sra-search 
	${SYS_LIBRARIES}
)

install ( TARGETS sra-search RUNTIME DESTINATION bin )

add_subdirectory ( test )
