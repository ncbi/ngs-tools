#set(CMAKE_VERBOSE_MAKEFILE on)

add_executable        ( compute-coverage compute-coverage.cpp )
add_dependencies      ( compute-coverage ncbi-vdb )
target_link_libraries ( compute-coverage
    tk-version
    ${CMAKE_STATIC_LIBRARY_PREFIX}kapp${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SYS_LIBRARIES}
)
links_and_install     ( compute-coverage )

add_executable        ( compute-coverage-ngs compute-coverage-ngs.cpp )
add_dependencies      ( compute-coverage-ngs ncbi-vdb )
target_link_libraries ( compute-coverage-ngs
    tk-version
    ${CMAKE_STATIC_LIBRARY_PREFIX}kapp${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SYS_LIBRARIES}
)
links_and_install     ( compute-coverage-ngs )

