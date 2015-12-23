cmake_minimum_required (VERSION 2.8)

# best to use at the top
project ( ngs-tools )

# this has to follow the project() command!
include ( build/cmake/env.cmake )  

include(ExternalProject)

if (WIN32)

ExternalProject_Add ( ngs 
	SOURCE_DIR ${NGS_ROOT}
	GIT_REPOSITORY https://github.com/ncbi/ngs.git
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
	BUILD_COMMAND msbuild ${NGS_ROOT}/ngs-sdk/win/ngs-sdk.sln /p:NGS_OUTDIR=${OUTDIR}/ngs-sdk/ /m /p:Platform=x64 /p:Configuration=Debug
	      COMMAND msbuild ${NGS_ROOT}/ngs-sdk/win/ngs-sdk.sln /p:NGS_OUTDIR=${OUTDIR}/ngs-sdk/ /m /p:Platform=x64 /p:Configuration=Release 
          COMMAND ant -f ${NGS_ROOT}/ngs-java -Dbuild.dir=${OUTDIR}/ngs-java jar
    INSTALL_COMMAND ""
)

ExternalProject_Add ( ncbi-vdb 
    DEPENDS ngs
	SOURCE_DIR ${VDB_ROOT}
	GIT_REPOSITORY https://github.com/ncbi/ncbi-vdb.git
    UPDATE_COMMAND git checkout VDB-2660-support-for-sra-search
    CONFIGURE_COMMAND ""
	BUILD_COMMAND msbuild ${VDB_ROOT}/build/MSVC/2010/ncbi-vdb.sln /p:NGS_OUTDIR=${OUTDIR}/ngs-sdk/ /p:VDB_OUTDIR=${OUTDIR}/ncbi-vdb/ /m /p:Platform=x64 /p:Configuration=Debug 
	      COMMAND msbuild ${VDB_ROOT}/build/MSVC/2010/ncbi-vdb.sln /p:NGS_OUTDIR=${OUTDIR}/ngs-sdk/ /p:VDB_OUTDIR=${OUTDIR}/ncbi-vdb/ /m /p:Platform=x64 /p:Configuration=Release
    INSTALL_COMMAND ""
)

else()

    ExternalProject_Add ( ngs 
        SOURCE_DIR ${NGS_ROOT}
        GIT_REPOSITORY https://github.com/ncbi/ngs.git
        CONFIGURE_COMMAND ${NGS_ROOT}/configure --prefix=${CMAKE_INSTALL_PREFIX} --build-prefix=${OUTDIR} ${CONFIGURE_FLAGS}
        BUILD_COMMAND make -C ${NGS_ROOT}/ngs-sdk COMMAND make -C ${NGS_ROOT}/ngs-java 
        INSTALL_COMMAND make -C ${NGS_ROOT}/ngs-sdk install COMMAND make -C ${NGS_ROOT}/ngs-java install
    )

    ExternalProject_Add ( ncbi-vdb 
        DEPENDS ngs
        SOURCE_DIR ${VDB_ROOT}
        GIT_REPOSITORY https://github.com/ncbi/ncbi-vdb.git
        UPDATE_COMMAND git checkout VDB-2660-support-for-sra-search
        CONFIGURE_COMMAND ${VDB_ROOT}/configure --prefix=${CMAKE_INSTALL_PREFIX} --build-prefix=${OUTDIR} ${CONFIGURE_FLAGS}
        BUILD_COMMAND make -C ${VDB_ROOT}
        INSTALL_COMMAND make -C ${VDB_ROOT} install
    )
    
endif()

add_subdirectory(sra-search)
add_subdirectory(frag-test)
add_subdirectory(frag-test-py)

set ( CPACK_PACKAGE_NAME ngs-tools )
set ( CPACK_PACKAGE_VERSION 0.1 )
set ( CPACK_PACKAGE_CONTACT "sra-tools@ncbi.nlm.nih.gov" )
set ( CPACK_PACKAGE_DESCRIPTION_SUMMARY "some tools that use NGS" )

set ( CPACK_DEBIAN_PACKAGE_ARCHITECTURE amd64 )

include ( CPack )
