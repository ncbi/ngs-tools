if( NOT WIN32 )

    macro( LoadDiffTest TestCase )
        add_test( NAME Test_LoadDiff_latf_sharq_${TestCase}
            COMMAND load-diff.sh ${SRATOOLS_BINDIR} ./temp
                "latf-load --quality PHRED_33 ./input/1.fastq -o ${TestCase}.latf" ${TestCase}.latf
                "sharq ./input/${TestCase}.fastq | ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/general-loader --schema sra/generic-fastq.vschema --target ${TestCase}.sharq" ${TestCase}.sharq
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )
    endmacro()

    LoadDiffTest( 1 )
    LoadDiffTest( 2 )

endif()