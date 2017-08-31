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

links_and_install ( sra-search )

add_subdirectory ( test )
