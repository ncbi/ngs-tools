add_executable ( general-loader
    general-loader.cpp
    protocol-parser.cpp
    database-loader.cpp
    main.cpp
)

include_directories ( ${CMAKE_SOURCE_DIR}/libs/general-writer/inc )

add_dependencies ( general-loader general-writer )

target_link_libraries ( general-loader
    general-writer
    ncbi-wvdb-static
    load
    kapp
    tk-version
	${SYS_LIBRARIES}
)

links_and_install ( general-loader )

add_subdirectory ( test )
