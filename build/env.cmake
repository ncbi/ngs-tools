#/////////////////////// Cache variables, may be overridden at config time:

if (UNIX)

    set ( PLATFORM x86_64 )

    if ( "${CMAKE_SYSTEM_NAME}" MATCHES "Darwin" )
        set ( OS mac )
        set ( COMPILER clang )
    else ()
        set ( OS linux )
        set ( COMPILER gcc )
    endif()

    # gmake is a single-configuration generator; we are either Debug or Release
    if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
        set ( BUILD dbg )
    else ()
        set ( BUILD rel )
    endif ()

    # by default, look for sister repositories sources side by side with ngs-tools, binaries under $HOME/ncbi-outdir
    set ( HOME "$ENV{HOME}" )

    set ( NGS_INCDIR  ${CMAKE_SOURCE_DIR}/../ngs/ngs-sdk                                              CACHE PATH "ngs include directory" )
    set ( NGS_LIBDIR  ${HOME}/ncbi-outdir/ngs-sdk/${OS}/${COMPILER}/${PLATFORM}/${BUILD}/lib          CACHE PATH "ngs library directory" )
    set ( VDB_INCDIR  ${CMAKE_SOURCE_DIR}/../ncbi-vdb/interfaces/                                     CACHE PATH "ncbi-vdb include directory" )
    set ( VDB_LIBDIR  ${HOME}/ncbi-outdir/ncbi-vdb/${OS}/${COMPILER}/${PLATFORM}/${BUILD}/lib         CACHE PATH "ncbi-vdb library directory" )
    set ( VDB_ILIBDIR ${HOME}/ncbi-outdir/ncbi-vdb/${OS}/${COMPILER}/${PLATFORM}/${BUILD}/ilib        CACHE PATH "ncbi-vdb internal library directory" )
    set ( SRATOOLS_BINDIR ${HOME}/ncbi-outdir/sra-tools/${OS}/${COMPILER}/${PLATFORM}/${BUILD}/bin    CACHE PATH "sra-tools executables directory" )

    set ( NGSTOOLS_OUTDIR ${HOME}/ncbi-outdir/ngs-tools/${OS}/${COMPILER}/${PLATFORM}/${BUILD} )

elseif (WIN32)

    set ( PLATFORM x64 )
    set ( OS win )
    set ( COMPILER "vc++" )
    if ( CMAKE_GENERATOR MATCHES "2010" )
        set ( PLATFORM_TOOLSET "v100" )
    elseif (CMAKE_GENERATOR MATCHES "2013" )
        set ( PLATFORM_TOOLSET "v120" )
    else()
        message( FATAL_ERROR "Unsupported generator for Windows: ${CMAKE_GENERATOR}." )
    endif()



    # by default, look for sister repositories sources side by side with ngs-tools, binaries under ../OUTDIR
    set ( NGS_INCDIR  ${CMAKE_SOURCE_DIR}/../ngs/ngs-sdk/                                                               CACHE PATH "ngs include directory" )
    set ( NGS_LIBDIR  ${CMAKE_SOURCE_DIR}/../OUTDIR/ngs-sdk/${OS}/${COMPILER}/${PLATFORM}/$(Configuration)/lib          CACHE PATH "ngs library directory" )
    set ( VDB_INCDIR  ${CMAKE_SOURCE_DIR}/../ncbi-vdb/interfaces/                                                       CACHE PATH "ncbi-vdb include directory" )
    set ( VDB_LIBDIR  ${CMAKE_SOURCE_DIR}/../OUTDIR/ncbi-vdb/${OS}/${COMPILER}/${PLATFORM}/$(Configuration)/lib         CACHE PATH "ncbi-vdb library directory" )
    set ( VDB_ILIBDIR ${CMAKE_SOURCE_DIR}/../OUTDIR/ncbi-vdb/${OS}/${COMPILER}/${PLATFORM}/$(Configuration)/ilib        CACHE PATH "ncbi-vdb internal library directory" )
    set ( SRATOOLS_BINDIR ${CMAKE_SOURCE_DIR}/../OUTDIR/sra-tools/${OS}/${COMPILER}/${PLATFORM}/$(Configuration)/bin    CACHE PATH "sra-tools executables directory" )

    set ( NGSTOOLS_OUTDIR ${CMAKE_SOURCE_DIR}/../OUTDIR/ngs-tools/${OS}/${COMPILER}/${PLATFORM}/$(Configuration) )

endif()

if ( NOT DEFINED CMAKE_RUNTIME_OUTPUT_DIRECTORY )
    set ( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${NGSTOOLS_OUTDIR}/bin )
endif()
if ( NOT DEFINED  CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
    set ( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${NGSTOOLS_OUTDIR}/ilib )
endif()

#/////////////////////////////////////////////////////////////////////////////////////////////

if (UNIX)

    if ( "${CMAKE_SYSTEM_NAME}" MATCHES "Darwin" )
        set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mmacosx-version-min=10.6 -stdlib=libstdc++" )
        # on Mac, we may need some gcc headers in addition to clang's
        include_directories ("${VDB_INCDIR}/cc/gcc/${PLATFORM}")
        include_directories ("${VDB_INCDIR}/cc/gcc")
    endif()

    include_directories ("${VDB_INCDIR}/os/unix")

    set ( SYS_LIBRARIES
            ${CMAKE_STATIC_LIBRARY_PREFIX}ncbi-ngs-c++${CMAKE_STATIC_LIBRARY_SUFFIX}
            ${CMAKE_STATIC_LIBRARY_PREFIX}ncbi-vdb-static${CMAKE_STATIC_LIBRARY_SUFFIX}
            ${CMAKE_STATIC_LIBRARY_PREFIX}ngs-c++${CMAKE_STATIC_LIBRARY_SUFFIX}
            pthread
            dl
    )

    if ( NOT DEFINED CMAKE_INSTALL_PREFIX)
        set ( CMAKE_INSTALL_PREFIX /usr/local/ )
    endif ()

    set ( CPACK_GENERATOR "RPM;DEB;TGZ;" )

elseif (WIN32)
    # Debug|Release is already in the NGSTOOLS_OUTDIR path; we do not want CMake to append /Debug|/Release to it
    SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG   ${NGSTOOLS_OUTDIR}/bin)
    SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${NGSTOOLS_OUTDIR}/bin)
    SET( CMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG   ${NGSTOOLS_OUTDIR}/ilib)
    SET( CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${NGSTOOLS_OUTDIR}/ilib)

    include_directories ("${NGS_INCDIR}/win")

    # use miltiple processors
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")

    # use Unicode
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /D_UNICODE /DUNICODE")

    # static run time libraries
    set ( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT" )
    set ( CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG}   /MTd" )

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

# look for dependencies

if (NOT EXISTS ${NGS_INCDIR})
    message( FATAL_ERROR "NGS includes are not found in ${NGS_INCDIR}." )
else()
    message( STATUS "Found NGS includes in ${NGS_INCDIR}. Looking for NGS libraries..." )

    if (UNIX)
        find_library ( NGS_LIBRARY ngs-c++ PATHS ${NGS_LIBDIR} NO_DEFAULT_PATH )
    else()
        if (CMAKE_CONFIGURATION_TYPES MATCHES "Debug")
            find_library ( NGS_LIBRARY libngs-bind-c++ PATHS ${NGS_LIBDIR_DEBUG} NO_DEFAULT_PATH )
        endif()
        if (NOT DEFINED NGS_LIBRARY AND CMAKE_CONFIGURATION_TYPES MATCHES "Release")
            find_library ( NGS_LIBRARY libngs-bind-c++ PATHS ${NGS_LIBDIR_RELEASE} NO_DEFAULT_PATH )
        endif()
    endif()

    if ( NGS_LIBRARY )
        get_filename_component(NGS_LIBRARY_DIR ${NGS_LIBRARY} PATH)
        message ( STATUS "Found NGS libraries in ${NGS_LIBDIR}" )
    else ()
        message( FATAL_ERROR "NGS libraries are not found in ${NGS_LIBDIR}." )
    endif()

    unset ( NGS_LIBRARY )
endif()

if (NOT EXISTS ${VDB_INCDIR})
    message( FATAL_ERROR "NCBI-VDB includes are not found in ${VDB_INCDIR}" )
else ()
    message( STATUS "Found NCBI-VDB includes in ${VDB_INCDIR}. Looking for NCBI-VDB libraries..." )
    if (UNIX)
        find_library ( VDB_LIBRARY ncbi-vdb PATHS ${VDB_LIBDIR} NO_DEFAULT_PATH )
    else()
        if ( CMAKE_CONFIGURATION_TYPES MATCHES "Debug" )
            find_library ( VDB_LIBRARY ncbi-vdb PATHS ${VDB_LIBDIR_DEBUG} NO_DEFAULT_PATH )
        endif()
        if ( CMAKE_CONFIGURATION_TYPES MATCHES "Release" )
            find_library ( VDB_LIBRARY ncbi-vdb PATHS ${VDB_LIBDIR_RELEASE} NO_DEFAULT_PATH )
        endif()
    endif()

    if ( VDB_LIBRARY )
        get_filename_component(VDB_LIBRARY_DIR ${VDB_LIBRARY} PATH)
        message ( STATUS "Found NGS libraries in ${VDB_LIBDIR}" )
    else()
        message( FATAL_ERROR "VDB libraries are not found in ${VDB_LIBDIR}." )
    endif()

    unset ( VDB_LIBRARY )
endif()

#//////////////////////////////////////////////////

include_directories ("${VDB_INCDIR}")
include_directories ("${VDB_INCDIR}/cc/${COMPILER}/${PLATFORM}")
include_directories ("${VDB_INCDIR}/cc/${COMPILER}")
include_directories ("${VDB_INCDIR}/os/${OS}")
include_directories ("${NGS_INCDIR}")

link_directories (  ${VDB_ILIBDIR} ${VDB_LIBDIR} ${NGS_LIBDIR} )

#/////////////////////////////////////////////////
# versioned names, symbolic links and installation for the tools

function ( links_and_install_subdir TARGET INST_SUBDIR)

    set_target_properties(${TARGET} PROPERTIES OUTPUT_NAME "${TARGET}.${VERSION}")

    set (TGTPATH ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${TARGET})

    ADD_CUSTOM_TARGET(link_${TARGET} ALL
                    COMMAND ${CMAKE_COMMAND}   -E create_symlink ${TARGET}.${VERSION} ${TGTPATH}.${MAJVERS}
                    COMMAND ${CMAKE_COMMAND}   -E create_symlink ${TARGET}.${MAJVERS} ${TGTPATH}
                    DEPENDS ${TARGET}
                    )

    install ( TARGETS ${TARGET} RUNTIME     DESTINATION bin/${INST_SUBDIR} )
    install ( FILES ${TGTPATH}.${MAJVERS}   DESTINATION bin/${INST_SUBDIR} )
    install ( FILES ${TGTPATH}              DESTINATION bin/${INST_SUBDIR} )

endfunction(links_and_install_subdir)

function ( links_and_install TARGET )
    links_and_install_subdir( ${TARGET} "" )
endfunction(links_and_install)

# testing
enable_testing()

