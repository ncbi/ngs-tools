cmake_minimum_required (VERSION 2.8)

find_package(Java REQUIRED)
include(UseJava)

include ( ${CMAKE_SOURCE_DIR}/build/cmake/env.cmake ) 

project ( ngs-fragtest )

add_jar ( ngs-fragtest
	SOURCES FragTest.java	     
	INCLUDE_JARS ${NGSJAR}
)

install_jar ( ngs-fragtest jar )
