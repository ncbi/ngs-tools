cmake_minimum_required (VERSION 2.8)

include ( ${CMAKE_SOURCE_DIR}/build/cmake/env.cmake ) 

project ( sra-search )

add_executable ( sra-search main.cpp vdb-search.cpp )

target_link_libraries ( sra-search 
	${SYS_LIBRARIES}
)

install ( TARGETS sra-search RUNTIME DESTINATION bin )

add_subdirectory ( test )
