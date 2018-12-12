include_directories ("..")

# command line tests
if ( NOT WIN32)
	set ( DUMP_REF_FASTA  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/dump-ref-fasta )

    # NB. Slow (run time > 1m) tests should names that start with "SlowTest_"; run without: "ctest -E SlowTest_", run only: "ctest -R SlowTest_"

	# help
	add_test ( NAME dump-ref-fasta-1.0-h WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
				${DUMP_REF_FASTA} 1.0-h --args "-h" )
	add_test ( NAME dump-ref-fasta-1.1-help WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
			   	${DUMP_REF_FASTA} 1.1-help --args "--help" )
	# version
    add_test ( NAME dump-ref-fasta-1.2-V WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND ./execandgrep.sh "${DUMP_REF_FASTA} -V" ${VERSION} )
    add_test ( NAME dump-ref-fasta-1.3-v WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND ./execandgrep.sh "${DUMP_REF_FASTA} -v" ${VERSION} )
    add_test ( NAME dump-ref-fasta-1.4-version WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND ./execandgrep.sh "${DUMP_REF_FASTA} --version" ${VERSION} )
	# regular run
	add_test ( NAME dump-ref-fasta-2.0 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh ${DUMP_REF_FASTA} 2.0 --args "SRR821492 | head -n 100" )
	# slice
    add_test ( NAME dump-ref-fasta-3.0 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh ${DUMP_REF_FASTA} 3.0 --args "SRR821492 NC_000002.11:1000-2000" )
	# local references only
    add_test ( NAME dump-ref-fasta-4.0 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh ${DUMP_REF_FASTA} 4.0 --args "--localref SRR821492" )
	# local references only, slice
    add_test ( NAME dump-ref-fasta-4.1 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh ${DUMP_REF_FASTA} 4.1 --args "-l SRR821492 NC_000002.11:1000-2000" )

endif ()
