add_subdirectory(compute-coverage)
add_subdirectory(ref-variation)
add_subdirectory(fastrq-dump)
add_subdirectory(sra-search)
add_subdirectory(sra_dnld_mgr)
add_subdirectory(dump_ref_fasta)

#add_subdirectory(compute-coverage) depends on VDB-3373 : will be added when VDB-3373 ncbi-vdb branch is merged to engineering

if ( NOT WIN32)
    add_subdirectory(general-loader)
    add_subdirectory(pileup-stats)
endif()

if (NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
   add_subdirectory(tax)
endif()

