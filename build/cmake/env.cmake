#/////////////////////// Cache variables, may be overridden at config time:

if (UNIX)
	set ( OLD_OUTDIR $ENV{HOME}/ncbi-outdir )
elseif (WIN32)
	set ( OLD_OUTDIR ${CMAKE_SOURCE_DIR}/../OUTDIR )
else ()
    message ( FATAL_ERROR "Unsupported OS" )
endif()
set ( CMAKE_OUTDIR ${CMAKE_BINARY_DIR} )

# by default, look for sister repositories sources side by side with ngs-tools
set ( NGS_ROOT  ${CMAKE_SOURCE_DIR}/../ngs            CACHE PATH "ngs source directory" )
set ( VDB_ROOT  ${CMAKE_SOURCE_DIR}/../ncbi-vdb       CACHE PATH "ncbi-vdb source directory" )
set ( NGS_SDK_BUILD_PREFIX ${OLD_OUTDIR}/ngs-sdk      CACHE PATH "ngs build directory" )
set ( VDB_BUILD_PREFIX ${OLD_OUTDIR}/ncbi-vdb         CACHE PATH "ncbi-vdb build directory" )
set ( NGS_TOOLS_OUTDIR_ENABLED OFF                    CACHE BOOL "enable copying build results into NGS_TOOLS_OUTDIR_PREFIX" )
set ( NGS_TOOLS_OUTDIR_PREFIX ${OLD_OUTDIR}/ngs-tools CACHE PATH "ngs-tools output directory" )

set ( GIT_BRANCH "master" )
set ( GIT_BRANCH_NGS ${GIT_BRANCH} CACHE STRING "git branch to use for ngs repository" )
set ( GIT_BRANCH_VDB ${GIT_BRANCH} CACHE STRING "git branch to use for ncbi-vdb repository" )

set ( PLATFORM x86_64 )
set ( WIN_PLATFORM x64 )

#/////////////////////////////////////////////////////////////////////////////////////////////

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

    set ( NGS_LIBDIR ${NGS_SDK_BUILD_PREFIX}/${OS}/${COMPILER}/${PLATFORM}/${BUILD}/lib )
    set ( VDB_LIBDIR ${VDB_BUILD_PREFIX}/${OS}/${COMPILER}/${PLATFORM}/${BUILD}/lib )
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

    # use miltiple processors
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")

    # use Unicode
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /D_UNICODE /DUNICODE")

    # static run time libraries
    set ( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT" )
    set ( CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG}   /MTd" )

    set ( NGS_LIBDIR ${NGS_SDK_BUILD_PREFIX}/${OS}/${PLATFORM_TOOLSET}/${WIN_PLATFORM}/$(Configuration)/lib )
    set ( VDB_LIBDIR ${VDB_BUILD_PREFIX}/${OS}/${PLATFORM_TOOLSET}/${WIN_PLATFORM}/$(Configuration)/lib )
    set ( VDB_ILIBDIR ${VDB_LIBDIR} )

    string(REPLACE "$(Configuration)" "Debug"   NGS_LIBDIR_DEBUG ${NGS_LIBDIR})
    string(REPLACE "$(Configuration)" "Release" NGS_LIBDIR_RELEASE ${NGS_LIBDIR})
    string(REPLACE "$(Configuration)" "Debug"   VDB_LIBDIR_DEBUG ${VDB_LIBDIR})
    string(REPLACE "$(Configuration)" "Release" VDB_LIBDIR_RELEASE ${VDB_LIBDIR})

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

# Java needs
find_package(Java REQUIRED)
find_program(ANT_EXECUTABLE ant PATHS $ENV{ANT_HOME} ENV PATH )
if ( NOT ANT_EXECUTABLE )
    message ( WARNING "Failed to locate 'ant' executable in PATH or ANT_HOME. Please specify path to 'ant' via ANT_EXECUTABLE" )
endif()
if ( NOT DEFINED ENV{JAVA_HOME} )
    message ( STATUS "Warning: JAVA_HOME is not set, 'ant' scripts may work incorrectly" )
endif ()
set ( NGSJAR "${OLD_OUTDIR}/ngs-java/jar/ngs-java.jar" )
set ( CMAKE_JAVA_COMPILE_FLAGS "-Xmaxerrs" "1" )

include ( functions )

#//////////////////////// External projects
# if project directory exists, do nothing, assume the user builds the project manually.
# otherwise, download and build the project

include(ExternalProject)

set ( EXTERNAL_PROJECTS "" )
set ( MSBUILD_OPTIONS  /m /tv:${TOOLS_VERSION} /p:Platform=${WIN_PLATFORM} )

if (NOT EXISTS ${NGS_ROOT})
    message( STATUS "NGS sources are not found. Will download and build them" )
    if (WIN32)
        ExternalProject_Add ( ngs
            SOURCE_DIR ${NGS_ROOT}
            GIT_REPOSITORY https://github.com/ncbi/ngs.git
            GIT_TAG ${GIT_BRANCH_NGS}
            UPDATE_COMMAND ""
            CONFIGURE_COMMAND ""
            BUILD_COMMAND msbuild ${MSBUILD_OPTIONS} ${NGS_ROOT}/ngs-sdk/win/${NGS_VSPROJ_SUBDIR}/ngs-sdk.sln /p:NGS_OUTDIR=${NGS_SDK_BUILD_PREFIX}/ /p:Configuration=Debug
                  COMMAND msbuild ${MSBUILD_OPTIONS} ${NGS_ROOT}/ngs-sdk/win/${NGS_VSPROJ_SUBDIR}/ngs-sdk.sln /p:NGS_OUTDIR=${NGS_SDK_BUILD_PREFIX}/ /p:Configuration=Release
                  COMMAND ${ANT_EXECUTABLE} -f ${NGS_ROOT}/ngs-java -Dbuild.dir=${OLD_OUTDIR}/ngs-java jar
          	INSTALL_COMMAND ""
        )
    else ()
        ExternalProject_Add ( ngs
            SOURCE_DIR ${NGS_ROOT}
            GIT_REPOSITORY https://github.com/ncbi/ngs.git
            GIT_TAG ${GIT_BRANCH_NGS}
            UPDATE_COMMAND ""
            CONFIGURE_COMMAND ${NGS_ROOT}/configure --prefix=${CMAKE_INSTALL_PREFIX} --build-prefix=${OLD_OUTDIR} ${CONFIGURE_FLAGS}
            BUILD_COMMAND make -C ${NGS_ROOT}/ngs-sdk COMMAND make -C ${NGS_ROOT}/ngs-java
            INSTALL_COMMAND ""
        )
    endif()
else()
    message( STATUS "Found NGS sources. Looking for NGS libraries..." )
    if (UNIX)
        find_library ( NGS_LIBRARY ngs-c++ PATHS ${NGS_LIBDIR} NO_DEFAULT_PATH )
    elseif (CMAKE_CONFIGURATION_TYPES MATCHES "Debug")
        find_library ( NGS_LIBRARY libngs-bind-c++ PATHS ${NGS_LIBDIR_DEBUG} NO_DEFAULT_PATH )
    else ()
        find_library ( NGS_LIBRARY libngs-bind-c++ PATHS ${NGS_LIBDIR_RELEASE} NO_DEFAULT_PATH )
    endif()
    if ( NGS_LIBRARY )
        get_filename_component(NGS_LIBRARY_DIR ${NGS_LIBRARY} PATH)
        message ( STATUS "Found NGS libraries at ${NGS_LIBRARY_DIR}" )
    else ()
        message ( STATUS "Warning: NGS libraries are not found at ${NGS_LIBDIR}. Will try to build them" )
    endif()
    unset ( NGS_LIBRARY CACHE )

    if (WIN32)
        ExternalProject_Add ( ngs
            SOURCE_DIR ${NGS_ROOT}
            CONFIGURE_COMMAND ""
            BUILD_COMMAND msbuild ${MSBUILD_OPTIONS} ${NGS_ROOT}/ngs-sdk/win/${NGS_VSPROJ_SUBDIR}/ngs-sdk.sln /p:NGS_OUTDIR=${NGS_SDK_BUILD_PREFIX}/ /p:Configuration=Debug
                  COMMAND msbuild ${MSBUILD_OPTIONS} ${NGS_ROOT}/ngs-sdk/win/${NGS_VSPROJ_SUBDIR}/ngs-sdk.sln /p:NGS_OUTDIR=${NGS_SDK_BUILD_PREFIX}/ /p:Configuration=Release
                  COMMAND ${ANT_EXECUTABLE} -f ${NGS_ROOT}/ngs-java -Dbuild.dir=${OLD_OUTDIR}/ngs-java jar
          	INSTALL_COMMAND ""
        )
    else ()
        ExternalProject_Add ( ngs
            SOURCE_DIR ${NGS_ROOT}
            CONFIGURE_COMMAND ""
            BUILD_COMMAND make -C ${NGS_ROOT}/ngs-sdk COMMAND make -C ${NGS_ROOT}/ngs-java
            INSTALL_COMMAND ""
        )
    endif()
endif()

if (NOT EXISTS ${VDB_ROOT})
    message( STATUS "NCBI-VDB sources are not found. Will download and build them." )
    if (WIN32)
        ExternalProject_Add ( ncbi-vdb
            DEPENDS ngs
            SOURCE_DIR ${VDB_ROOT}
            GIT_REPOSITORY https://github.com/ncbi/ncbi-vdb.git
            GIT_TAG ${GIT_BRANCH_VDB}
            UPDATE_COMMAND ""
            CONFIGURE_COMMAND ""
            BUILD_COMMAND msbuild ${MSBUILD_OPTIONS} ${VDB_ROOT}/build/MSVC/${VDB_VSPROJ_SUBDIR}/ncbi-vdb.vcxproj   /p:VDB_OUTDIR=${OLD_OUTDIR}/ /p:Configuration=Debug
                  COMMAND msbuild ${MSBUILD_OPTIONS} ${VDB_ROOT}/build/MSVC/${VDB_VSPROJ_SUBDIR}/ncbi-vdb.vcxproj   /p:VDB_OUTDIR=${OLD_OUTDIR}/ /p:Configuration=Release
                  COMMAND msbuild ${MSBUILD_OPTIONS} ${VDB_ROOT}/build/MSVC/${VDB_VSPROJ_SUBDIR}/bz2.vcxproj        /p:VDB_OUTDIR=${OLD_OUTDIR}/ /p:Configuration=Debug
                  COMMAND msbuild ${MSBUILD_OPTIONS} ${VDB_ROOT}/build/MSVC/${VDB_VSPROJ_SUBDIR}/bz2.vcxproj        /p:VDB_OUTDIR=${OLD_OUTDIR}/ /p:Configuration=Release
                  COMMAND msbuild ${MSBUILD_OPTIONS} ${VDB_ROOT}/build/MSVC/${VDB_VSPROJ_SUBDIR}/zlib.vcxproj       /p:VDB_OUTDIR=${OLD_OUTDIR}/ /p:Configuration=Debug
                  COMMAND msbuild ${MSBUILD_OPTIONS} ${VDB_ROOT}/build/MSVC/${VDB_VSPROJ_SUBDIR}/zlib.vcxproj       /p:VDB_OUTDIR=${OLD_OUTDIR}/ /p:Configuration=Release
                  COMMAND msbuild ${MSBUILD_OPTIONS} ${VDB_ROOT}/build/MSVC/${VDB_VSPROJ_SUBDIR}/ngs-c++.vcxproj    /p:VDB_OUTDIR=${OLD_OUTDIR}/ /p:Configuration=Debug
                  COMMAND msbuild ${MSBUILD_OPTIONS} ${VDB_ROOT}/build/MSVC/${VDB_VSPROJ_SUBDIR}/ngs-c++.vcxproj    /p:VDB_OUTDIR=${OLD_OUTDIR}/ /p:Configuration=Release
                  COMMAND msbuild ${MSBUILD_OPTIONS} ${VDB_ROOT}/build/MSVC/${VDB_VSPROJ_SUBDIR}/ktst.vcxproj       /p:VDB_OUTDIR=${OLD_OUTDIR}/ /p:Configuration=Debug
                  COMMAND msbuild ${MSBUILD_OPTIONS} ${VDB_ROOT}/build/MSVC/${VDB_VSPROJ_SUBDIR}/ktst.vcxproj       /p:VDB_OUTDIR=${OLD_OUTDIR}/ /p:Configuration=Release
                  COMMAND msbuild ${MSBUILD_OPTIONS} ${VDB_ROOT}/build/MSVC/${VDB_VSPROJ_SUBDIR}/kapp.vcxproj       /p:VDB_OUTDIR=${OLD_OUTDIR}/ /p:Configuration=Debug
                  COMMAND msbuild ${MSBUILD_OPTIONS} ${VDB_ROOT}/build/MSVC/${VDB_VSPROJ_SUBDIR}/kapp.vcxproj       /p:VDB_OUTDIR=${OLD_OUTDIR}/ /p:Configuration=Release
                  COMMAND msbuild ${MSBUILD_OPTIONS} ${VDB_ROOT}/build/MSVC/${VDB_VSPROJ_SUBDIR}/load.vcxproj       /p:VDB_OUTDIR=${OLD_OUTDIR}/ /p:Configuration=Debug
                  COMMAND msbuild ${MSBUILD_OPTIONS} ${VDB_ROOT}/build/MSVC/${VDB_VSPROJ_SUBDIR}/load.vcxproj       /p:VDB_OUTDIR=${OLD_OUTDIR}/ /p:Configuration=Release
          	INSTALL_COMMAND ""
        )
    else()
        ExternalProject_Add ( ncbi-vdb
            DEPENDS ngs
            SOURCE_DIR ${VDB_ROOT}
            GIT_REPOSITORY https://github.com/ncbi/ncbi-vdb.git
            GIT_TAG ${GIT_BRANCH_VDB}
            UPDATE_COMMAND ""
            CONFIGURE_COMMAND ${VDB_ROOT}/configure --prefix=${CMAKE_INSTALL_PREFIX} --build-prefix=${OLD_OUTDIR} ${CONFIGURE_FLAGS}
            BUILD_COMMAND make -C ${VDB_ROOT}
            INSTALL_COMMAND ""
        )
    endif()
else ()
    message( STATUS "Found NCBI-VDB sources. Looking for NCBI-VDB libraries..." )
    if (UNIX)
        find_library ( VDB_LIBRARY ncbi-vdb PATHS ${VDB_LIBDIR} NO_DEFAULT_PATH )
    elseif (CMAKE_CONFIGURATION_TYPES MATCHES "Debug")
        find_library ( VDB_LIBRARY ncbi-vdb PATHS ${VDB_LIBDIR_DEBUG} NO_DEFAULT_PATH )
    else ()
        find_library ( VDB_LIBRARY ncbi-vdb PATHS ${VDB_LIBDIR_RELEASE} NO_DEFAULT_PATH )
    endif()

    if ( VDB_LIBRARY )
        get_filename_component(VDB_LIBRARY_DIR ${VDB_LIBRARY} PATH)
        message ( STATUS "Found NGS libraries at ${VDB_LIBRARY_DIR}" )
    else ()
        message ( STATUS "Warning: VDB libraries are not found at ${VDB_LIBDIR}. Will try to build them" )
    endif()
    unset ( VDB_LIBRARY CACHE )

    if (WIN32)
        ExternalProject_Add ( ncbi-vdb
            DEPENDS ngs
            SOURCE_DIR ${VDB_ROOT}
            CONFIGURE_COMMAND ""
            BUILD_COMMAND msbuild ${MSBUILD_OPTIONS} ${VDB_ROOT}/build/MSVC/${VDB_VSPROJ_SUBDIR}/ncbi-vdb.vcxproj   /p:VDB_OUTDIR=${OLD_OUTDIR}/ /p:Configuration=Debug
                  COMMAND msbuild ${MSBUILD_OPTIONS} ${VDB_ROOT}/build/MSVC/${VDB_VSPROJ_SUBDIR}/ncbi-vdb.vcxproj   /p:VDB_OUTDIR=${OLD_OUTDIR}/ /p:Configuration=Release
                  COMMAND msbuild ${MSBUILD_OPTIONS} ${VDB_ROOT}/build/MSVC/${VDB_VSPROJ_SUBDIR}/bz2.vcxproj        /p:VDB_OUTDIR=${OLD_OUTDIR}/ /p:Configuration=Debug
                  COMMAND msbuild ${MSBUILD_OPTIONS} ${VDB_ROOT}/build/MSVC/${VDB_VSPROJ_SUBDIR}/bz2.vcxproj        /p:VDB_OUTDIR=${OLD_OUTDIR}/ /p:Configuration=Release
                  COMMAND msbuild ${MSBUILD_OPTIONS} ${VDB_ROOT}/build/MSVC/${VDB_VSPROJ_SUBDIR}/zlib.vcxproj       /p:VDB_OUTDIR=${OLD_OUTDIR}/ /p:Configuration=Debug
                  COMMAND msbuild ${MSBUILD_OPTIONS} ${VDB_ROOT}/build/MSVC/${VDB_VSPROJ_SUBDIR}/zlib.vcxproj       /p:VDB_OUTDIR=${OLD_OUTDIR}/ /p:Configuration=Release
                  COMMAND msbuild ${MSBUILD_OPTIONS} ${VDB_ROOT}/build/MSVC/${VDB_VSPROJ_SUBDIR}/ngs-c++.vcxproj    /p:VDB_OUTDIR=${OLD_OUTDIR}/ /p:Configuration=Debug
                  COMMAND msbuild ${MSBUILD_OPTIONS} ${VDB_ROOT}/build/MSVC/${VDB_VSPROJ_SUBDIR}/ngs-c++.vcxproj    /p:VDB_OUTDIR=${OLD_OUTDIR}/ /p:Configuration=Release
                  COMMAND msbuild ${MSBUILD_OPTIONS} ${VDB_ROOT}/build/MSVC/${VDB_VSPROJ_SUBDIR}/ktst.vcxproj       /p:VDB_OUTDIR=${OLD_OUTDIR}/ /p:Configuration=Debug
                  COMMAND msbuild ${MSBUILD_OPTIONS} ${VDB_ROOT}/build/MSVC/${VDB_VSPROJ_SUBDIR}/ktst.vcxproj       /p:VDB_OUTDIR=${OLD_OUTDIR}/ /p:Configuration=Release
          	INSTALL_COMMAND ""
        )
    else()
        ExternalProject_Add ( ncbi-vdb
            DEPENDS ngs
            SOURCE_DIR ${VDB_ROOT}
            CONFIGURE_COMMAND ""
            BUILD_COMMAND make -C ${VDB_ROOT}
            INSTALL_COMMAND ""
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

# testing
enable_testing()

