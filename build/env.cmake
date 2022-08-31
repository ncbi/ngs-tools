#/////////////////////// Cache variables, may be overridden at config time:

# by default, look for sister repositories sources side by side with ngs-tools, binaries under $OUTDIR if set, otherwise $HOME/ncbi-outdir
if (NOT DEFINED OUTDIR)
    if (DEFINED ENV{OUTDIR})
        set ( OUTDIR "$ENV{OUTDIR}" )
    else ()
        set ( OUTDIR "$ENV{HOME}/ncbi-outdir" )
    endif ()
endif()

if (UNIX)

#    set ( PLATFORM x86_64 )

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

    set ( VDB_INCDIR  ${CMAKE_SOURCE_DIR}/../ncbi-vdb/interfaces/                           CACHE PATH "ncbi-vdb include directory" )
    set ( VDB_LIBDIR  ${OUTDIR}/ncbi-vdb/${OS}/${COMPILER}/${PLATFORM}/${BUILD}/lib         CACHE PATH "ncbi-vdb library directory" )
    set ( VDB_ILIBDIR ${OUTDIR}/ncbi-vdb/${OS}/${COMPILER}/${PLATFORM}/${BUILD}/ilib        CACHE PATH "ncbi-vdb internal library directory" )
    set ( SRATOOLS_BINDIR ${OUTDIR}/sra-tools/${OS}/${COMPILER}/${PLATFORM}/${BUILD}/bin    CACHE PATH "sra-tools executables directory" )
    set ( NGSTOOLS_OUTDIR ${OUTDIR}/ngs-tools/${OS}/${COMPILER}/${PLATFORM}/${BUILD}        CACHE PATH "ngs-tools output directory" )

elseif (WIN32)

    set ( PLATFORM "x64" CACHE STRING "Windows Platform (x64 or Win32)" )
    set ( OS win )
    set ( COMPILER "vc++" )
    if ( CMAKE_GENERATOR MATCHES "2010" )
        set ( PLATFORM_TOOLSET "v100" )
    elseif (CMAKE_GENERATOR MATCHES "2013" )
        set ( PLATFORM_TOOLSET "v120" )
    elseif (CMAKE_GENERATOR MATCHES "2017" )
        set ( PLATFORM_TOOLSET "v141" )
    else()
        message( FATAL_ERROR "Unsupported generator for Windows: ${CMAKE_GENERATOR}." )
    endif()

    # by default, look for sister repositories sources side by side with ngs-tools, binaries under ../OUTDIR
    set ( VDB_INCDIR  ${CMAKE_SOURCE_DIR}/../ncbi-vdb/interfaces/                                           CACHE PATH "ncbi-vdb include directory" )
    set ( VDB_LIBDIR  ${OUTDIR}/ncbi-vdb/${OS}/${PLATFORM_TOOLSET}/${PLATFORM}/$(Configuration)/lib         CACHE PATH "ncbi-vdb library directory" )
    set ( VDB_ILIBDIR ${OUTDIR}/ncbi-vdb/${OS}/${PLATFORM_TOOLSET}/${PLATFORM}/$(Configuration)/ilib        CACHE PATH "ncbi-vdb internal library directory" )
    set ( SRATOOLS_BINDIR ${OUTDIR}/sra-tools/${OS}/${PLATFORM_TOOLSET}/${PLATFORM}/$(Configuration)/bin    CACHE PATH "sra-tools executables directory" )
    set ( NGSTOOLS_OUTDIR ${OUTDIR}/ngs-tools/${OS}/${COMPILER}/${PLATFORM}/${BUILD}                        CACHE PATH "ngs-tools output directory")

endif()

#/////////////////////////////////////////////////////////////////////////////////////////////

if (UNIX)

    # default executables and libaries output directories
    if ( NOT DEFINED CMAKE_RUNTIME_OUTPUT_DIRECTORY )
        set ( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${NGSTOOLS_OUTDIR}/bin )
    endif()
    if ( NOT DEFINED  CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
        set ( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${NGSTOOLS_OUTDIR}/ilib )
    endif()

    if ( "${CMAKE_SYSTEM_NAME}" MATCHES "Darwin" )
        set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mmacosx-version-min=10.10 -stdlib=libc++" )
        # on Mac, we may need some gcc headers in addition to clang's
        include_directories ("${VDB_INCDIR}/cc/gcc/${PLATFORM}")
        include_directories ("${VDB_INCDIR}/cc/gcc")
    endif()

    include_directories ("${VDB_INCDIR}/os/unix")

    set ( SYS_LIBRARIES
            ${CMAKE_STATIC_LIBRARY_PREFIX}ncbi-vdb-static${CMAKE_STATIC_LIBRARY_SUFFIX}
            pthread
            dl
    )

    set ( SYS_WLIBRARIES
            ${CMAKE_STATIC_LIBRARY_PREFIX}ncbi-wvdb-static${CMAKE_STATIC_LIBRARY_SUFFIX}
            pthread
            dl
    )

    if ( NOT DEFINED CMAKE_INSTALL_PREFIX)
        set ( CMAKE_INSTALL_PREFIX /usr/local/ )
    endif ()

    set ( CPACK_GENERATOR "RPM;DEB;TGZ;" )

elseif (WIN32)
    if ( NOT DEFINED CMAKE_RUNTIME_OUTPUT_DIRECTORY )
        set ( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${NGSTOOLS_OUTDIR}/${OS}/${PLATFORM_TOOLSET}/${PLATFORM} )
    endif()
    if ( NOT DEFINED  CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
        set ( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${NGSTOOLS_OUTDIR}/${OS}/${PLATFORM_TOOLSET}/${PLATFORM} )
    endif()
    # Work configuration name into the NGSTOOLS_OUTDIR path; we do not want CMake to add /Debug|/Release at the end
    SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG   ${NGSTOOLS_OUTDIR}/${OS}/${PLATFORM_TOOLSET}/${PLATFORM}/Debug/bin)
    SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${NGSTOOLS_OUTDIR}/${OS}/${PLATFORM_TOOLSET}/${PLATFORM}/Release/bin)
    SET( CMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG   ${NGSTOOLS_OUTDIR}/${OS}/${PLATFORM_TOOLSET}/${PLATFORM}/Debug/ilib)
    SET( CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${NGSTOOLS_OUTDIR}/${OS}/${PLATFORM_TOOLSET}/${PLATFORM}/Release/ilib)

    # use miltiple processors
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")

    # use Unicode
    add_definitions(-DUNICODE -D_UNICODE)

    # static run time libraries
    set ( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT" )
    set ( CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG}   /MTd" )
    set ( CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /MT" )
    set ( CMAKE_C_FLAGS_DEBUG   "${CMAKE_C_FLAGS_DEBUG}   /MTd" )

    string(REPLACE "$(Configuration)" "Debug"   VDB_LIBDIR_DEBUG ${VDB_LIBDIR})
    string(REPLACE "$(Configuration)" "Release" VDB_LIBDIR_RELEASE ${VDB_LIBDIR})

    set ( SYS_LIBRARIES
        ${CMAKE_STATIC_LIBRARY_PREFIX}bz2${CMAKE_STATIC_LIBRARY_SUFFIX}
        ${CMAKE_STATIC_LIBRARY_PREFIX}zlib${CMAKE_STATIC_LIBRARY_SUFFIX}
        ${CMAKE_STATIC_LIBRARY_PREFIX}ncbi-vdb${CMAKE_STATIC_LIBRARY_SUFFIX}
        ws2_32
        Crypt32
    )

    set ( CPACK_GENERATOR "ZIP" )

else()
    message ( FATAL_ERROR "Unsupported OS" )
endif()

# Java needs
find_package(Java)
if ( Java_FOUND )
    find_program(ANT_EXECUTABLE ant PATHS $ENV{ANT_HOME} ENV PATH )
    if ( NOT ANT_EXECUTABLE )
        message ( WARNING "Failed to locate 'ant' executable in PATH or ANT_HOME. Please specify path to 'ant' via ANT_EXECUTABLE" )
    endif()
    if ( NOT DEFINED ENV{JAVA_HOME} )
        message ( STATUS "Warning: JAVA_HOME is not set, 'ant' scripts may work incorrectly" )
    endif ()
    set ( CMAKE_JAVA_COMPILE_FLAGS "-Xmaxerrs" "1" )
endif()

# look for dependencies
if (NOT EXISTS ${VDB_INCDIR})
    message( FATAL_ERROR "NCBI-VDB includes are not found in ${VDB_INCDIR}" )
else ()
    message( STATUS "Found NCBI-VDB includes in ${VDB_INCDIR}. Looking for NCBI-VDB libraries..." )
    if (UNIX)
        find_library ( VDB_LIBRARY ncbi-vdb PATHS ${VDB_LIBDIR} NO_DEFAULT_PATH )
    else()

		# on Windows, require both debug and release libraries
        if (CMAKE_CONFIGURATION_TYPES MATCHES ".*Debug.*")
            find_library ( VDB_LIBRARY ncbi-vdb PATHS ${VDB_LIBDIR_DEBUG} NO_DEFAULT_PATH )
			if ( VDB_LIBRARY )
				get_filename_component(VDB_LIBRARY ${VDB_LIBRARY} PATH)
				message ( STATUS "Found Debug NCBI-VDB libraries in ${VDB_LIBRARY}" )
			else ()
				message( FATAL_ERROR "NCBI-VDB libraries are not found in ${VDB_LIBDIR_DEBUG}." )
			endif()
        endif()

        if (CMAKE_CONFIGURATION_TYPES MATCHES ".*Release.*")
            find_library ( VDB_LIBRARY ncbi-vdb PATHS ${VDB_LIBDIR_RELEASE} NO_DEFAULT_PATH )
			if ( VDB_LIBRARY )
				get_filename_component(VDB_LIBRARY ${VDB_LIBRARY} PATH)
				message ( STATUS "Found Release NCBI-VDB libraries in ${VDB_LIBDIR_RELEASE}" )
			else ()
				message( FATAL_ERROR "NCBI-VDB libraries are not found in ${VDB_LIBDIR_RELEASE}." )
			endif()
        endif()

    endif()

    unset ( VDB_LIBRARY )
endif()

#//////////////////////////////////////////////////

include_directories ("${VDB_INCDIR}")
include_directories ("${VDB_INCDIR}/cc/${COMPILER}/${PLATFORM}")
include_directories ("${VDB_INCDIR}/cc/${COMPILER}")
include_directories ("${VDB_INCDIR}/os/${OS}")

link_directories (  ${VDB_ILIBDIR} ${VDB_LIBDIR} )

#/////////////////////////////////////////////////
# versioned names, symbolic links and installation for the tools

function ( links_and_install_subdir TARGET INST_SUBDIR)

    if (WIN32)

        install ( TARGETS ${TARGET} RUNTIME DESTINATION bin )

    else()

        # on Unix, version the binary and add 2 symbolic links
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

    endif()

endfunction(links_and_install_subdir)

function ( links_and_install TARGET )
    links_and_install_subdir( ${TARGET} "" )
endfunction(links_and_install)

# testing
enable_testing()

