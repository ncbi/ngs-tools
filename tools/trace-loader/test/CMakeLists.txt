include_directories ( ".." ${CMAKE_SOURCE_DIR}/libs )

if ( NOT WIN32)
    ## Hola mothra, hello bathra
    find_program( TRA_LD trace-loader HINTS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} )
    if ( NOT TRA_LD )
        message( SEND_ERROR "trace-loader is not found, exiting with error" )
    endif ()
    message( STATUS "trace-loader found in ${TRA_LD}" )

    find_program( GEN_LD general-loader HINTS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} )
    if ( NOT GEN_LD )
        message( SEND_ERROR "general-loader is not found, exiting with error" )
    endif ()
    message( STATUS "general-loader found in ${GEN_LD}" )

    find_program( VDB_DM vdb-dump HINTS ${SRATOOLS_BINDIR} )
    if ( NOT VDB_DM )
        message( SEND_ERROR "vdb-dump is not found, exiting with error" )
    endif ()
    message( STATUS "vdb-dump found in ${VDB_DM}" )

    add_custom_command (
            OUTPUT test/run-test.log
            COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/run-test.sh ${CMAKE_CURRENT_SOURCE_DIR} ${TRA_LD} ${GEN_LD} ${VDB_DM} ${VDB_INCDIR} ${CMAKE_CURRENT_SOURCE_DIR}/../rules
            COMMENT "Holla bolla"
            )

    add_custom_target(
            RTS ALL
            DEPENDS test/run-test.log
            )
endif ()
