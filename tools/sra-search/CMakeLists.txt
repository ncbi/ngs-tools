add_executable ( sra-search
    main.cpp
    searchblock.cpp
    fragmentmatchiterator.cpp
    blobmatchiterator.cpp
    referencematchiterator.cpp
    vdb-search.cpp
)

include_directories ( ${CMAKE_SOURCE_DIR}/libs )

add_dependencies ( sra-search ngs-vdb )

target_link_libraries ( sra-search
    ngs-vdb
	${SYS_LIBRARIES}
)

ncbi_copy_exec_to_old_outdir ( sra-search )
install ( TARGETS sra-search RUNTIME DESTINATION bin )

add_subdirectory ( test )
