### ref-variation

add_executable ( ref-variation common.cpp helper.cpp ref-variation.cpp  )

ncbi_copy_exec_to_old_outdir ( ref-variation )

# add_dependencies ( ref-variation ngs-vdb )

target_link_libraries ( ref-variation
    ngs-vdb
    ${CMAKE_STATIC_LIBRARY_PREFIX}kapp${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${CMAKE_STATIC_LIBRARY_PREFIX}load${CMAKE_STATIC_LIBRARY_SUFFIX}
	${SYS_LIBRARIES}
)

### var-expand

add_executable ( var-expand common.cpp helper.cpp var-expand.cpp )

ncbi_copy_exec_to_old_outdir ( var-expand )

target_link_libraries ( var-expand
    ${CMAKE_STATIC_LIBRARY_PREFIX}kapp${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${CMAKE_STATIC_LIBRARY_PREFIX}load${CMAKE_STATIC_LIBRARY_SUFFIX}
	${SYS_LIBRARIES}
)

### common

include_directories ( ${CMAKE_SOURCE_DIR}/libs )

# linker flags
if (WIN32)
    set ( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /ENTRY:wmainCRTStartup")
endif()

install ( TARGETS ref-variation var-expand RUNTIME DESTINATION bin )

add_subdirectory ( test )

