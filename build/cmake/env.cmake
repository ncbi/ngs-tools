#/////////////////////// Cache variables, may be overridden at config time:

set ( OUTDIR ${CMAKE_BINARY_DIR} )

# by default, look for sister repositories sources side by side with ngs-tools
set ( NGS_ROOT  ${CMAKE_SOURCE_DIR}/../ngs            CACHE PATH "ngs source directory" )
set ( VDB_ROOT  ${CMAKE_SOURCE_DIR}/../ncbi-vdb        CACHE PATH "ncbi-vdb source directory" )

set ( INSTALL_DIR ${CMAKE_SOURCE_DIR}/../INSTALL    CACHE PATH "installation directory" )

set ( GIT_BRANCH "master" )
set ( GIT_BRANCH_NGS ${GIT_BRANCH} CACHE STRING "git branch to use for ngs repository" )
set ( GIT_BRANCH_VDB ${GIT_BRANCH} CACHE STRING "git branch to use for ncbi-vdb repository" )

#/////////////////////////////////////////////////////////////////////////////////////////////

set ( PLATFORM x86_64 )

if (UNIX)

    if ( "${CMAKE_SYSTEM_NAME}" MATCHES "Darwin" )
        set ( OS mac )
        set ( COMPILER clang )
        set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mmacosx-version-min=10.6 -stdlib=libstdc++" )
        # on Mac, we may need some gcc headers in addition to clang's
        include_directories ("${VDB_ROOT}/interfaces/cc/gcc/${PLATFORM}")
        include_directories ("${VDB_ROOT}/interfaces/cc/gcc")
    else ()
        set ( OS linux )
        set ( COMPILER gcc )
    endif()

    include_directories ("${VDB_ROOT}/interfaces/os/unix")
    
    # gmake is a single-configuration generator; we are either Debug or Release
    if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
        set ( BUILD dbg )
        set ( CONFIGURE_FLAGS "--with-debug" )
    else ()
        set ( BUILD rel )
        set ( CONFIGURE_FLAGS "" )
    endif ()
    
    set ( NGS_LIBDIR ${OUTDIR}/ngs-sdk/${OS}/${COMPILER}/${PLATFORM}/${BUILD}/lib )
    set ( VDB_LIBDIR ${OUTDIR}/ncbi-vdb/${OS}/${COMPILER}/${PLATFORM}/${BUILD}/lib )
    set ( VDB_ILIBDIR ${VDB_LIBDIR}/../ilib/ )

    set ( SYS_LIBRARIES 
            ${CMAKE_STATIC_LIBRARY_PREFIX}ncbi-ngs-c++${CMAKE_STATIC_LIBRARY_SUFFIX}
            ${CMAKE_STATIC_LIBRARY_PREFIX}ncbi-vdb-static${CMAKE_STATIC_LIBRARY_SUFFIX}
            ${CMAKE_STATIC_LIBRARY_PREFIX}ngs-c++${CMAKE_STATIC_LIBRARY_SUFFIX}
            pthread 
            dl 
    )
    
    if (!CMAKE_INSTALL_PREFIX)
        set ( CMAKE_INSTALL_PREFIX /usr/local/ )  
    endif ()

    set ( CPACK_GENERATOR "RPM;DEB;TGZ;" )
    
elseif (WIN32)

    set ( OS win )
    set ( COMPILER vc++ )

    if ( CMAKE_GENERATOR MATCHES "2010" )
        set ( NGS_VSPROJ_SUBDIR "vs2010" ) 
        set ( VDB_VSPROJ_SUBDIR "2010" ) 
        set ( PLATFORM_TOOLSET "v100" ) 
        set ( TOOLS_VERSION "4.0" )
    elseif (CMAKE_GENERATOR MATCHES "2013" )
        set ( NGS_VSPROJ_SUBDIR "vs2013" ) 
        set ( VDB_VSPROJ_SUBDIR "2013" ) 
        set ( PLATFORM_TOOLSET "v120" ) 
        set ( TOOLS_VERSION "12.0" )
    else()
        message( FATAL_ERROR "Unsupported generator for Windows: ${CMAKE_GENERATOR}." )        
    endif()

    include_directories ("${NGS_ROOT}/ngs-sdk/win")
    
    set ( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT" )
    set ( CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG}   /MTd" )

    set ( NGS_LIBDIR ${OUTDIR}/ngs-sdk/${OS}/${PLATFORM_TOOLSET}/$(Platform)/$(Configuration)/lib )
    set ( VDB_LIBDIR ${OUTDIR}/ncbi-vdb/${OS}/${PLATFORM_TOOLSET}/$(Platform)/$(Configuration)/lib )
    set ( VDB_ILIBDIR ${VDB_LIBDIR} )

    set ( SYS_LIBRARIES 
        ${CMAKE_STATIC_LIBRARY_PREFIX}bz2${CMAKE_STATIC_LIBRARY_SUFFIX}
        ${CMAKE_STATIC_LIBRARY_PREFIX}zlib${CMAKE_STATIC_LIBRARY_SUFFIX}
        ${CMAKE_STATIC_LIBRARY_PREFIX}ngs-c++${CMAKE_STATIC_LIBRARY_SUFFIX}
        ${CMAKE_STATIC_LIBRARY_PREFIX}ncbi-vdb${CMAKE_STATIC_LIBRARY_SUFFIX}
        libngs-bind-c++${CMAKE_STATIC_LIBRARY_SUFFIX}
        libngs-disp${CMAKE_STATIC_LIBRARY_SUFFIX}
        ws2_32
    )

    set ( CPACK_GENERATOR "ZIP" )
    
else()

    message ( FATAL_ERROR "Unsupported OS" )
    
endif()

#//////////////////////// External projects
# if project directory exists, do nothing, assume the user builds the project manually.
# otherwise, download and build the project

include(ExternalProject)

set ( EXTERNAL_PROJECTS "" )

if (NOT EXISTS ${NGS_ROOT})
    if (WIN32)
        ExternalProject_Add ( ngs 
            SOURCE_DIR ${NGS_ROOT}
            GIT_REPOSITORY https://github.com/ncbi/ngs.git
            GIT_TAG ${GIT_BRANCH_NGS}
            UPDATE_COMMAND ""
            CONFIGURE_COMMAND ""
            BUILD_COMMAND msbuild /m ${NGS_ROOT}/ngs-sdk/win/${NGS_VSPROJ_SUBDIR}/ngs-sdk.sln /tv:${TOOLS_VERSION} /p:NGS_OUTDIR=${OUTDIR}/ngs-sdk/ /m /p:Platform=x64 /p:Configuration=Debug
                  COMMAND msbuild /m ${NGS_ROOT}/ngs-sdk/win/${NGS_VSPROJ_SUBDIR}/ngs-sdk.sln /tv:${TOOLS_VERSION} /p:NGS_OUTDIR=${OUTDIR}/ngs-sdk/ /m /p:Platform=x64 /p:Configuration=Release 
                  COMMAND ant -f ${NGS_ROOT}/ngs-java -Dbuild.dir=${OUTDIR}/ngs-java jar
            INSTALL_COMMAND ""
        )
    else ()
        ExternalProject_Add ( ngs 
            SOURCE_DIR ${NGS_ROOT}
            GIT_REPOSITORY https://github.com/ncbi/ngs.git
            GIT_TAG ${GIT_BRANCH_NGS}
            UPDATE_COMMAND ""
            CONFIGURE_COMMAND ${NGS_ROOT}/configure --prefix=${CMAKE_INSTALL_PREFIX} --build-prefix=${OUTDIR} ${CONFIGURE_FLAGS}
            BUILD_COMMAND make -C ${NGS_ROOT}/ngs-sdk COMMAND make -C ${NGS_ROOT}/ngs-java 
            INSTALL_COMMAND make -C ${NGS_ROOT}/ngs-sdk install COMMAND make -C ${NGS_ROOT}/ngs-java install
        )
    endif()
else()
    if (WIN32)
        ExternalProject_Add ( ngs 
            SOURCE_DIR ${NGS_ROOT}
            CONFIGURE_COMMAND ""
            BUILD_COMMAND msbuild /m ${NGS_ROOT}/ngs-sdk/win/${NGS_VSPROJ_SUBDIR}/ngs-sdk.sln /tv:${TOOLS_VERSION} /p:NGS_OUTDIR=${OUTDIR}/ngs-sdk/ /m /p:Platform=x64 /p:Configuration=Debug
                  COMMAND msbuild /m ${NGS_ROOT}/ngs-sdk/win/${NGS_VSPROJ_SUBDIR}/ngs-sdk.sln /tv:${TOOLS_VERSION} /p:NGS_OUTDIR=${OUTDIR}/ngs-sdk/ /m /p:Platform=x64 /p:Configuration=Release 
                  COMMAND ant -f ${NGS_ROOT}/ngs-java -Dbuild.dir=${OUTDIR}/ngs-java jar
            INSTALL_COMMAND ""
        )
    else ()
        ExternalProject_Add ( ngs 
            SOURCE_DIR ${NGS_ROOT}
            CONFIGURE_COMMAND ${NGS_ROOT}/configure --prefix=${CMAKE_INSTALL_PREFIX} --build-prefix=${OUTDIR} ${CONFIGURE_FLAGS}
            BUILD_COMMAND make -C ${NGS_ROOT}/ngs-sdk COMMAND make -C ${NGS_ROOT}/ngs-java 
            INSTALL_COMMAND make -C ${NGS_ROOT}/ngs-sdk install COMMAND make -C ${NGS_ROOT}/ngs-java install
        )
    endif()
endif()

if (NOT EXISTS ${VDB_ROOT})
    if (WIN32)
        ExternalProject_Add ( ncbi-vdb 
            DEPENDS ngs
            SOURCE_DIR ${VDB_ROOT}
            GIT_REPOSITORY https://github.com/ncbi/ncbi-vdb.git
            GIT_TAG ${GIT_BRANCH_VDB}
            UPDATE_COMMAND ""
            CONFIGURE_COMMAND ""
            BUILD_COMMAND msbuild /m ${VDB_ROOT}/build/MSVC/${VDB_VSPROJ_SUBDIR}/ncbi-vdb.vcxproj /tv:${TOOLS_VERSION} /p:NGS_OUTDIR=${OUTDIR}/ngs-sdk/ /p:VDB_OUTDIR=${OUTDIR}/ncbi-vdb/ /m /p:Platform=x64 /p:Configuration=Debug 
                  COMMAND msbuild /m ${VDB_ROOT}/build/MSVC/${VDB_VSPROJ_SUBDIR}/ncbi-vdb.vcxproj /tv:${TOOLS_VERSION} /p:NGS_OUTDIR=${OUTDIR}/ngs-sdk/ /p:VDB_OUTDIR=${OUTDIR}/ncbi-vdb/ /m /p:Platform=x64 /p:Configuration=Release
            INSTALL_COMMAND ""
        )
    else()
        ExternalProject_Add ( ncbi-vdb 
            DEPENDS ngs
            SOURCE_DIR ${VDB_ROOT}
            GIT_REPOSITORY https://github.com/ncbi/ncbi-vdb.git
            GIT_TAG ${GIT_BRANCH_VDB}
            UPDATE_COMMAND ""
            CONFIGURE_COMMAND ${VDB_ROOT}/configure --prefix=${CMAKE_INSTALL_PREFIX} --build-prefix=${OUTDIR} ${CONFIGURE_FLAGS}
            BUILD_COMMAND make -C ${VDB_ROOT}
            INSTALL_COMMAND make -C ${VDB_ROOT} install
        )
    endif()
else ()
    if (WIN32)
        ExternalProject_Add ( ncbi-vdb 
            DEPENDS ngs
            SOURCE_DIR ${VDB_ROOT}
            CONFIGURE_COMMAND ""
            BUILD_COMMAND msbuild /m ${VDB_ROOT}/build/MSVC/${VDB_VSPROJ_SUBDIR}/ncbi-vdb.vcxproj /tv:${TOOLS_VERSION} /p:NGS_OUTDIR=${OUTDIR}/ngs-sdk/ /p:VDB_OUTDIR=${OUTDIR}/ncbi-vdb/ /m /p:Platform=x64 /p:Configuration=Debug 
                  COMMAND msbuild /m ${VDB_ROOT}/build/MSVC/${VDB_VSPROJ_SUBDIR}/ncbi-vdb.vcxproj /tv:${TOOLS_VERSION} /p:NGS_OUTDIR=${OUTDIR}/ngs-sdk/ /p:VDB_OUTDIR=${OUTDIR}/ncbi-vdb/ /m /p:Platform=x64 /p:Configuration=Release
            INSTALL_COMMAND ""
        )
    else()
        ExternalProject_Add ( ncbi-vdb 
            DEPENDS ngs
            SOURCE_DIR ${VDB_ROOT}
            CONFIGURE_COMMAND ${VDB_ROOT}/configure --prefix=${CMAKE_INSTALL_PREFIX} --build-prefix=${OUTDIR} ${CONFIGURE_FLAGS}
            BUILD_COMMAND make -C ${VDB_ROOT}
            INSTALL_COMMAND make -C ${VDB_ROOT} install
        )
    endif()
endif()

#//////////////////////////////////////////////////

include_directories ("${VDB_ROOT}/interfaces")
include_directories ("${VDB_ROOT}/interfaces/cc/${COMPILER}/${PLATFORM}")
include_directories ("${VDB_ROOT}/interfaces/cc/${COMPILER}")
include_directories ("${VDB_ROOT}/interfaces/os/${OS}")
include_directories ("${NGS_ROOT}/ngs-sdk")

link_directories (  ${VDB_ILIBDIR} ${VDB_LIBDIR} ${NGS_LIBDIR} )

# Java needs
set ( NGSJAR "${OUTDIR}/ngs-java/jar/ngs-java.jar" )
set ( CMAKE_JAVA_COMPILE_FLAGS "-Xmaxerrs" "1" )

# testing
enable_testing()

