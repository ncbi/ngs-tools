set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /SUBSYSTEM:WINDOWS /ENTRY:wmainCRTStartup")
add_executable ( fastrq-dump fastrq-dump.cpp filter.cpp formatter.cpp operation.cpp settings.cpp )

target_link_libraries ( fastrq-dump
    tk-version
    ${CMAKE_STATIC_LIBRARY_PREFIX}kapp${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SYS_LIBRARIES}
)

links_and_install ( fastrq-dump )