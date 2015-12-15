cmake_minimum_required (VERSION 2.8)

include ( ${CMAKE_SOURCE_DIR}/build/cmake/env.cmake ) 

project ( sra-search )

add_executable ( sra-search main.cpp vdb-search.cpp )

target_link_libraries ( sra-search 
	${VDB_LIBDIR}/${CMAKE_STATIC_LIBRARY_PREFIX}ncbi-ngs-c++${CMAKE_STATIC_LIBRARY_SUFFIX}
	${VDB_LIBDIR}/${CMAKE_STATIC_LIBRARY_PREFIX}ncbi-vdb${CMAKE_STATIC_LIBRARY_SUFFIX}
	${SYS_LIBRARIES}
)

install ( TARGETS sra-search RUNTIME DESTINATION bin )
 

