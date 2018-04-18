#set(CMAKE_VERBOSE_MAKEFILE on)
if (WIN32)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /SUBSYSTEM:WINDOWS /ENTRY:wmainCRTStartup")
endif ()

add_executable        ( compute-coverage compute-coverage.cpp )
target_link_libraries ( compute-coverage
    tk-version
    ${CMAKE_STATIC_LIBRARY_PREFIX}kapp${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SYS_LIBRARIES}
)
links_and_install     ( compute-coverage )

add_executable        ( compute-coverage-ngs compute-coverage-ngs.cpp )
target_link_libraries ( compute-coverage-ngs
    tk-version
    ${CMAKE_STATIC_LIBRARY_PREFIX}kapp${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SYS_LIBRARIES}
)
links_and_install     ( compute-coverage-ngs )

