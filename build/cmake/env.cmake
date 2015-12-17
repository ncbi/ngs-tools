set ( PLATFORM x86_64 )

if (UNIX)
	# TODO: support gmake on Mac
	
	set ( OS linux )
	set ( COMPILER gcc )

	set ( OUTDIR /home/ncbi/devel/OUTDIR )
	set ( NGS_ROOT /home/ncbi/devel/ngs )
	set ( VDB_ROOT /home/ncbi/devel/ncbi-vdb )

	# gmake is a single-configuration generator; we are either Debug or Release
	if (CMAKE_BUILD_TYPE STREQUAL "Debug")
		set ( BUILD dbg )
	else ()
		set ( BUILD rel )
	endif ()
	# select Debug/Release versions of our "3d party" libraries
	set ( NGS_LIBDIR ${OUTDIR}/ngs-sdk/${OS}/${COMPILER}/${PLATFORM}/${BUILD}/lib/ )
	set ( VDB_LIBDIR ${OUTDIR}/ncbi-vdb/${OS}/${COMPILER}/${PLATFORM}/${BUILD}/lib/ )
	set ( VDB_ILIBDIR ${OUTDIR}/ncbi-vdb/${OS}/${COMPILER}/${PLATFORM}/${BUILD}/ilib/ )

	set ( SYS_LIBRARIES 
			${NGS_LIBDIR}/${CMAKE_STATIC_LIBRARY_PREFIX}ngs-c++${CMAKE_STATIC_LIBRARY_SUFFIX}
			pthread 
			dl 
	)

	include_directories ("${VDB_ROOT}/interfaces/os/unix")

	set ( CMAKE_INSTALL_PREFIX /usr/local/ )  

	set ( CPACK_GENERATOR "RPM;DEB;TGZ;" )
	
elseif (WIN32)
	set ( OS win )
	set ( COMPILER vc++ )

	set ( OUTDIR C:/Users/NCBI/github/OUTDIR )
	set ( NGS_ROOT C:/Users/NCBI/github/ngs )
	set ( VDB_ROOT C:/Users/NCBI/github/ncbi-vdb )

	# switch between Debug/Release versions of our "3d party" libraries
	set ( BUILD $<CONFIG> )
	
	# TODO: adjust external builds to be consistent with Unix, then switch to PLATFORM and COMPILER here;
	# NGS_LIBDIR and VDB_LIBDIR then can be defined outside of this if
	set ( WIN_PLATFORM x64 )
	set ( WIN_COMPILER cl )
	set ( NGS_LIBDIR ${OUTDIR}/${OS}/${WIN_COMPILER}/${WIN_PLATFORM}/${BUILD}/lib/ )
	set ( VDB_LIBDIR ${OUTDIR}/${OS}/${WIN_COMPILER}/${WIN_PLATFORM}/${BUILD}/lib/ )
	set ( VDB_ILIBDIR ${VDB_LIBDIR} )

	set ( SYS_LIBRARIES 
		${VDB_LIBDIR}/${CMAKE_STATIC_LIBRARY_PREFIX}bz2${CMAKE_STATIC_LIBRARY_SUFFIX}
		${VDB_LIBDIR}/${CMAKE_STATIC_LIBRARY_PREFIX}zlib${CMAKE_STATIC_LIBRARY_SUFFIX}
		${NGS_LIBDIR}/libngs-bind-c++${CMAKE_STATIC_LIBRARY_SUFFIX}
		${NGS_LIBDIR}/libngs-disp${CMAKE_STATIC_LIBRARY_SUFFIX}
		ws2_32
	)

	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")

	set ( CMAKE_INSTALL_PREFIX "C:/Program Files/ngs-tools" )
	  
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

# Java needs
set ( NGSJAR "${OUTDIR}/ngs-java/jar/ngs-java.jar" )
set ( CMAKE_JAVA_COMPILE_FLAGS "-Xmaxerrs" "1" )
