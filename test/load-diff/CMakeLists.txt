if( NOT WIN32 )

    add_test( NAME Test_LoadDiff_latf_sharq
        COMMAND load-diff.sh ${SRATOOLS_BINDIR} ${VDB_INCDIR} ./temp
            "latf-load --quality PHRED_33 ./input/1.fastq -o 1.latf" latf
            "sharq ./input/1.fastq | ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/general-loader --include ${VDB_INCDIR} --schema sra/generic-fastq.vschema --target 1.sharq --log-level=info --xml-log=log" sharq
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

endif()