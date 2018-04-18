set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /SUBSYSTEM:WINDOWS /ENTRY:wmainCRTStartup")
add_executable ( dump-ref-fasta DumpReferenceFASTA.cpp)

target_link_libraries ( dump-ref-fasta
    tk-version
    ${CMAKE_STATIC_LIBRARY_PREFIX}kapp${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SYS_LIBRARIES}
)

links_and_install (dump-ref-fasta)
