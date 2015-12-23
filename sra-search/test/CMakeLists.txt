cmake_minimum_required (VERSION 2.8)

project ( test-sra-search )

include_directories ("..")

add_executable ( test-sra-search test-sra-search.cpp )
add_dependencies ( test-sra-search sra-search ngs ncbi-vdb ) 

target_link_libraries ( test-sra-search 
	${SYS_LIBRARIES}
    ${CMAKE_STATIC_LIBRARY_PREFIX}ktst${CMAKE_STATIC_LIBRARY_SUFFIX}
)

# white box test
add_test ( SraSearchWhiteBox test-sra-search )

# command line tests
if ( NOT WIN32)

	# bad args
	add_test ( NAME SraSearch-1.0-h 				WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh ${CMAKE_CURRENT_BINARY_DIR}/../sra-search 1.0-h  					--args "-h" ) 
	add_test ( NAME SraSearch-1.1-help 				WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh ${CMAKE_CURRENT_BINARY_DIR}/../sra-search 1.1-help  				--args "--help" ) 
	add_test ( NAME SraSearch-2.0-no-args 			WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh ${CMAKE_CURRENT_BINARY_DIR}/../sra-search 2.0-no-args 			--rc 1 ) 
	add_test ( NAME SraSearch-2.1-no-run 			WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh ${CMAKE_CURRENT_BINARY_DIR}/../sra-search 2.1-no-run 				--rc 1 --args AGCT ) 
	add_test ( NAME SraSearch-2.2-missing-algorithm WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh ${CMAKE_CURRENT_BINARY_DIR}/../sra-search 2.2-missing-algorithm 	--rc 1 --args "AGCT SRR000001 -a" )
	 
endif ()
