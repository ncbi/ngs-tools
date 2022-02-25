if( NOT WIN32 )

    # specify the location of schema files in a local .kfg file, to be used by the tests here as needed
    add_test(NAME Setup COMMAND bash -c "echo 'vdb/schema/paths = \"${VDB_INCDIR}\"' > tmp.kfg" WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )

    macro( LoadDiffTest TestCase )
        add_test( NAME Test_LoadDiff_latf_sharq_${TestCase}
            COMMAND
                ${CMAKE_COMMAND} -E env NCBI_SETTINGS=/
                ${CMAKE_COMMAND} -E env VDB_CONFIG=${CMAKE_CURRENT_SOURCE_DIR}
                ./load-diff.sh ${SRATOOLS_BINDIR} ./temp
                "latf-load --quality PHRED_33 ./input/${TestCase} -o ${TestCase}.latf" ${TestCase}.latf
                "sharq ./input/${TestCase} | ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/general-loader --schema sra/generic-fastq.vschema --target ${TestCase}.sharq  -L debug" ${TestCase}.sharq
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )
    endmacro()

    macro( LoadDiffTestPair TestCase File1 File2 latf-params sharq-params)
        add_test( NAME Test_LoadDiff_latf_sharq_${TestCase}
            COMMAND
                ${CMAKE_COMMAND} -E env NCBI_SETTINGS=/
                ${CMAKE_COMMAND} -E env VDB_CONFIG=${CMAKE_CURRENT_SOURCE_DIR}
                ./load-diff.sh ${SRATOOLS_BINDIR} ./temp
                "latf-load ${latf-params}  -1 ./input/${File1} -2 ./input/${File2} -o ${TestCase}.latf" ${TestCase}.latf
                "sharq ${sharq-params} ./input/${File1} ./input/${File2} | ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/general-loader --schema sra/generic-fastq.vschema --target ${TestCase}.sharq" ${TestCase}.sharq
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )
    endmacro()

    LoadDiffTest( 1.fastq )
    LoadDiffTest( 2.fastq )
    LoadDiffTest( BRA_AJ155.fastq )
    LoadDiffTestPair( SRR12060554_IR SRR12060554_I1.fastq SRR12060554_R2.fastq "--quality PHRED_33" "--quality 33" )
    LoadDiffTestPair( stLFR2_split_read stLFR2_split_read.1.fq stLFR2_split_read.2.fq "--quality PHRED_33" "--quality 33" )
    LoadDiffTestPair( skin_S1_L001 skin_S1_L001_I1_001.fastq skin_S1_L001_R2_001.fastq "--quality PHRED_33" "--quality 33" )
    LoadDiffTestPair( S1_3 S1_3_1.fq S1_3_2.fq "--quality PHRED_64" "--quality 64" )
    LoadDiffTest( glucose.RNA-seq.single.MinION.fastq )

endif()