include_directories ( ${CMAKE_SOURCE_DIR} )

add_library ( ngs-vdb
                VdbReadCollection.cpp
                FragmentBlobIterator.cpp
                FragmentBlob.cpp
                NGS-VDB.cpp )

add_subdirectory ( test )
