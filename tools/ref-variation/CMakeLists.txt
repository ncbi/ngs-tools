### ref-variation

add_executable ( ref-variation common.cpp helper.cpp ref-variation.cpp  )

# add_dependencies ( ref-variation ngs-vdb )

target_link_libraries ( ref-variation
    tk-version
    ngs-vdb
    ${CMAKE_STATIC_LIBRARY_PREFIX}kapp${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${CMAKE_STATIC_LIBRARY_PREFIX}load${CMAKE_STATIC_LIBRARY_SUFFIX}
	${SYS_LIBRARIES}
)

links_and_install ( ref-variation )

### var-expand

add_executable ( var-expand common.cpp helper.cpp var-expand.cpp )

target_link_libraries ( var-expand
    tk-version
    ${CMAKE_STATIC_LIBRARY_PREFIX}kapp${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${CMAKE_STATIC_LIBRARY_PREFIX}load${CMAKE_STATIC_LIBRARY_SUFFIX}
	${SYS_LIBRARIES}
)

links_and_install ( var-expand )

### common

include_directories ( ${CMAKE_SOURCE_DIR}/libs )

# linker flags
if (WIN32)
    set ( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /ENTRY:wmainCRTStartup")
endif()

add_subdirectory ( test )

