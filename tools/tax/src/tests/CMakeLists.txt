include_directories ( ${CMAKE_SOURCE_DIR} )
include_directories ( ${CMAKE_CURRENT_SOURCE_DIR}/.. )

add_executable ( hash           hash )
add_executable ( reader_test    reader_test ${CMAKE_CURRENT_SOURCE_DIR}/../reader.cpp )
add_executable ( seq_transform  seq_transform )

target_link_libraries ( hash ${SYS_LIBRARIES} )
target_link_libraries ( reader_test ${SYS_LIBRARIES} )
target_link_libraries ( seq_transform ${SYS_LIBRARIES} )

add_test ( NAME hash COMMAND hash )
add_test ( NAME reader_test COMMAND reader_test WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/.. )
add_test ( NAME seq_transform COMMAND seq_transform )
