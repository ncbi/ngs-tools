include_directories ( ${CMAKE_SOURCE_DIR}/libs/general-writer/inc )

# libgeneral-writer
add_library ( general-writer
                general-writer.cpp
                utf8-like-int-codec.c
             )

# gw-dumper

add_executable( gw-dumper gw-dumper.cpp )

add_dependencies ( gw-dumper general-writer )

target_link_libraries ( gw-dumper
    general-writer
	${SYS_LIBRARIES}
)
add_subdirectory ( test )
