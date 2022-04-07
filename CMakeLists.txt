cmake_minimum_required (VERSION 2.8.12)

message(CMAKE_VERSION = ${CMAKE_VERSION})

set ( CMAKE_OSX_DEPLOYMENT_TARGET 10.10 )

set ( CMAKE_CXX_EXTENSIONS OFF )
if (CMAKE_VERSION VERSION_LESS "3.1")
   set (CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")
   message(CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS})
else ()
   set (CMAKE_CXX_STANDARD 11)
   message(CMAKE_CXX_STANDARD = ${CMAKE_CXX_STANDARD})
endif ()


set(CMAKE_VERBOSE_MAKEFILE ON)

# best to use at the top
project ( ngs-tools )
if( PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR )
    message(FATAL_ERROR "In-source builds are not allowed")
endif()


# this has to follow the project() command!
set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    "${CMAKE_CURRENT_SOURCE_DIR}/build"
)
include ( env NO_POLICY_SCOPE )

add_subdirectory(libs)
add_subdirectory(tools)
add_subdirectory(test)

set ( CPACK_PACKAGE_NAME ngs-tools )
set ( CPACK_PACKAGE_VERSION 0.1 )
set ( CPACK_PACKAGE_CONTACT "sra-tools@ncbi.nlm.nih.gov" )
set ( CPACK_PACKAGE_DESCRIPTION_SUMMARY "some tools that use NGS" )

set ( CPACK_DEBIAN_PACKAGE_ARCHITECTURE amd64 )

include ( CPack )
