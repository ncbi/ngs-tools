cmake_minimum_required (VERSION 2.8)

# best to use at the top
project ( ngs-tools )

# this has to follow the project() command!
include ( build/cmake/env.cmake )  

add_subdirectory(sra-search)
add_subdirectory(tools/fastrq-dump)
add_subdirectory(sra_dnld_mgr)

set ( CPACK_PACKAGE_NAME ngs-tools )
set ( CPACK_PACKAGE_VERSION 0.1 )
set ( CPACK_PACKAGE_CONTACT "sra-tools@ncbi.nlm.nih.gov" )
set ( CPACK_PACKAGE_DESCRIPTION_SUMMARY "some tools that use NGS" )

set ( CPACK_DEBIAN_PACKAGE_ARCHITECTURE amd64 )

include ( CPack )
