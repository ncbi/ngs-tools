if( NOT WIN32 )

    # specify the location of schema files in a local .kfg file, to be used by the tests here as needed
    add_test(NAME Setup COMMAND bash -c "echo 'vdb/schema/paths = \"${VDB_INCDIR}\"' > tmp.kfg" WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )

    macro( SharqTest TestCase )
        add_test( NAME Test_LoadDiff_latf_sharq_${TestCase}
            COMMAND
                ${CMAKE_COMMAND} -E env NCBI_SETTINGS=/
                ${CMAKE_COMMAND} -E env VDB_CONFIG=${CMAKE_CURRENT_SOURCE_DIR}
                ./load-dump-diff.sh ${SRATOOLS_BINDIR} ./temp
                "sharq ./input/${TestCase} | ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/general-loader --schema sra/generic-fastq.vschema --target ${TestCase}  " ${TestCase}
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )
    endmacro()

    SharqTest( glucose.RNA-seq.single.MinION.fastq )

endif()