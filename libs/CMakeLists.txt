add_subdirectory(tk-version)
add_subdirectory(ngs-vdb)

if ( NOT WIN32)
    add_subdirectory(general-writer)
endif()
