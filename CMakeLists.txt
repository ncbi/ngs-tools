cmake_minimum_required (VERSION 2.8)

project ( ngs-tools )

add_subdirectory(sra-search)
add_subdirectory(frag-test)
add_subdirectory(frag-test-py)

include ( build/cmake/env.cmake )  

set ( CPACK_PACKAGE_NAME ngs-tools )
set ( CPACK_PACKAGE_VERSION 0.1 )
set ( CPACK_PACKAGE_CONTACT "sra-tools@ncbi.nlm.nih.gov" )
set ( CPACK_PACKAGE_DESCRIPTION_SUMMARY "some tools that use NGS" )

set ( CPACK_DEBIAN_PACKAGE_ARCHITECTURE amd64 )

include ( CPack )
