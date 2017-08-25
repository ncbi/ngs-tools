include_directories ( ${CMAKE_SOURCE_DIR}/libs/general-writer/inc )

add_library ( general-writer
                general-writer.cpp
                utf8-like-int-codec.c
             )

add_subdirectory ( test )
