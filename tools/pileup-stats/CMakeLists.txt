add_executable ( pileup-stats
    pileup-stats
)

include_directories ( ${CMAKE_SOURCE_DIR}/libs/general-writer/inc )

add_dependencies ( pileup-stats general-writer )

target_link_libraries ( pileup-stats
    general-writer
    ncbi-wvdb-static
    kapp
    tk-version
	${SYS_LIBRARIES}
)

install ( TARGETS pileup-stats RUNTIME DESTINATION bin )

add_subdirectory ( test )
