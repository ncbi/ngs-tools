set ( PLATFORM x86_64 )

set ( NGS_ROOT "$ENV{HOME}/ngs" )
set ( VDB_ROOT "$ENV{HOME}/ncbi-vdb" )
set ( INSTALL_DIR "$ENV{HOME}/install" )

# by default, look for sister repositories sources side by side with ngs-tools
if (NOT DEFINED NGS_ROOT)
	set ( NGS_ROOT  ${CMAKE_SOURCE_DIR}/../ngs )
endif()

if (NOT DEFINED VDB_ROOT)
	set ( VDB_ROOT  ${CMAKE_SOURCE_DIR}/../ncbi-vdb )
endif()

if (UNIX)

    set ( OUTDIR ${CMAKE_BINARY_DIR}/.. )
    
	# TODO: support gmake on Mac
	
	set ( OS linux )
	set ( COMPILER gcc )

	# gmake is a single-configuration generator; we are either Debug or Release
	if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
		set ( BUILD dbg )
        set ( CONFIGURE_FLAGS "--with-debug" )
	else ()
		set ( BUILD rel )
        set ( CONFIGURE_FLAGS "" )
	endif ()
    
	# By default, look for our "3d party" libraries side by side with our binary directory
	set ( NGS_LIBDIR ${NGS_ROOT}/../OUTDIR/ngs-sdk/${OS}/${COMPILER}/${PLATFORM}/${BUILD}/lib )
	set ( VDB_LIBDIR ${VDB_ROOT}/../OUTDIR/ncbi-vdb/${OS}/${COMPILER}/${PLATFORM}/${BUILD}/lib )
	set ( VDB_ILIBDIR ${VDB_LIBDIR}/../ilib/ )

	set ( SYS_LIBRARIES 
            ${CMAKE_STATIC_LIBRARY_PREFIX}ncbi-ngs-c++${CMAKE_STATIC_LIBRARY_SUFFIX}
            ${CMAKE_STATIC_LIBRARY_PREFIX}ncbi-vdb${CMAKE_STATIC_LIBRARY_SUFFIX}
			${CMAKE_STATIC_LIBRARY_PREFIX}ngs-c++${CMAKE_STATIC_LIBRARY_SUFFIX}
			pthread 
			dl 
	)
	include_directories ("${VDB_ROOT}/interfaces/os/unix")

    if (!CMAKE_INSTALL_PREFIX)
    	set ( CMAKE_INSTALL_PREFIX /usr/local/ )  
    endif ()

	set ( CPACK_GENERATOR "RPM;DEB;TGZ;" )
	
elseif (WIN32)

    set ( OUTDIR ${CMAKE_BINARY_DIR} )

	set ( OS win )
	set ( COMPILER vc++ )

	# By default, look for our "3d party" libraries side by side with our binary directory
	set ( NGS_LIBDIR ${OUTDIR}/../OUTDIR/${OS}/cl/$(Platform)/$(Configuration)/lib )
	set ( VDB_LIBDIR ${OUTDIR}/../OUTDIR/${OS}/cl/$(Platform)/$(Configuration)/lib )
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

	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
	set(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG}   /MTd")

    include_directories ("${NGS_ROOT}/ngs-sdk/win")
    
    if (!CMAKE_INSTALL_PREFIX)
        set ( CMAKE_INSTALL_PREFIX "C:/Program Files/ngs-tools" )
    endif ()
	  
	set ( CPACK_GENERATOR "ZIP" )
	
else()
	# TODO: support XCode on Mac

	message ( FATAL_ERROR "Unsupported OS" )
endif()


include_directories ("${VDB_ROOT}/interfaces")
include_directories ("${VDB_ROOT}/interfaces/cc/${COMPILER}/${PLATFORM}")
include_directories ("${VDB_ROOT}/interfaces/cc/${COMPILER}")
include_directories ("${VDB_ROOT}/interfaces/os/${OS}")
include_directories ("${NGS_ROOT}/ngs-sdk")

link_directories (  ${VDB_ILIBDIR} ${VDB_LIBDIR} ${NGS_LIBDIR} )

# Java needs
set ( NGSJAR "${NGS_ROOT}/../OUTDIR/ngs-java/jar/ngs-java.jar" )
set ( CMAKE_JAVA_COMPILE_FLAGS "-Xmaxerrs" "1" )

# testing
enable_testing()