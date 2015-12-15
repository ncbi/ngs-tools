cmake_minimum_required (VERSION 2.8)

include ( ${CMAKE_SOURCE_DIR}/build/cmake/env.cmake ) 

project ( other-app )

add_executable ( other-app main.cpp )

target_link_libraries ( other-app 
	${VDB_LIBDIR}/${CMAKE_STATIC_LIBRARY_PREFIX}ncbi-ngs-c++${CMAKE_STATIC_LIBRARY_SUFFIX}
	${VDB_LIBDIR}/${CMAKE_STATIC_LIBRARY_PREFIX}ncbi-vdb${CMAKE_STATIC_LIBRARY_SUFFIX}
	${SYS_LIBRARIES}
)

install ( TARGETS other-app RUNTIME DESTINATION bin )
 
