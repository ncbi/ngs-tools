# command line tests
if ( NOT WIN32)
	add_test ( NAME RunAll WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND run_all.sh )
endif ()
