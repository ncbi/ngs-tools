include_directories ( ${CMAKE_SOURCE_DIR}/libs/general-writer/inc )

add_executable ( test-general-writer test-general-writer.cpp )
add_dependencies ( test-general-writer general-writer )

target_link_libraries ( test-general-writer
    general-writer
    ${CMAKE_STATIC_LIBRARY_PREFIX}ktst${CMAKE_STATIC_LIBRARY_SUFFIX}
	${SYS_LIBRARIES}
)

add_test ( NAME GeneralWriterTest
           COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test-general-writer
           WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/libs/general-writer/src/test )
