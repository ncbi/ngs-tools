include_directories ("..")

add_executable ( test-sra-search test-sra-search.cpp ../vdb-search.cpp ../searchblock.cpp)
add_dependencies ( test-sra-search sra-search ngs-vdb )

target_link_libraries ( test-sra-search
    ngs-vdb
	${SYS_LIBRARIES}
    ${CMAKE_STATIC_LIBRARY_PREFIX}ktst${CMAKE_STATIC_LIBRARY_SUFFIX}
)

# white box test
add_test ( SraSearchWhiteBox test-sra-search )

# command line tests
if ( NOT WIN32)
	set ( SRA_SEARCH  ${CMAKE_CURRENT_BINARY_DIR}/../sra-search )

	# help
	add_test ( NAME SraSearch-1.0-h WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
				${SRA_SEARCH} 1.0-h --args "-h" )
	add_test ( NAME SraSearch-1.1-help WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
			   	${SRA_SEARCH} 1.1-help --args "--help" )

	# bad args, bad arguments
	add_test ( NAME SraSearch-2.0-no-args WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
				${SRA_SEARCH} 2.0-no-args --rc 1 )
	add_test ( NAME SraSearch-2.1-no-run WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
				${SRA_SEARCH} 2.1-no-run --rc 1 --args AGCT )
	add_test ( NAME SraSearch-2.2-missing-algorithm WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
				${SRA_SEARCH} 2.2-missing-algorithm --rc 1 --args "AGCT SRR000001 -a" )
	add_test ( NAME SraSearch-2.3-bad-algorithm WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
				${SRA_SEARCH} 2.3-bad-algorithm --rc 1 --args "--algorithm bad AGCT SRR000001" )

	# single run
	add_test ( NAME SraSearch-3.05-single-run-default-WGS WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
				${SRA_SEARCH} 3.05-single-run-default-WGS --args "AGCTAGCTA ALAI01" )

	# query expressions, bad arguments
	add_test ( NAME SraSearch-5.1-run-not-nucstrstr-with-expression WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
				${SRA_SEARCH} 5.1-run-not-nucstrstr-with-expression --args "AAAAAAACCCCCCCAAAAAAACCCCCCC\\|\\|AGCTAGCTAGCT --algorithm FgrepStandard --expression  SRR000001" --rc 1 )

	# imperfect match, bad arguments
	add_test ( NAME SraSearch-6.0.1-imperfect-match-no-arg WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
				${SRA_SEARCH} 6.0.1-imperfect-match-no-arg --args "AGCTAGCTAGCT SRR000001 -S" --rc 1 )
	add_test ( NAME SraSearch-6.0.2-imperfect-match-bad-arg WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
				${SRA_SEARCH} 6.0.2-imperfect-match-bad-arg --args "AGCTAGCTAGCT SRR000001 -S qq" --rc 1  )
	add_test ( NAME SraSearch-6.0.3-imperfect-match-default WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
	         	${SRA_SEARCH} 6.0.3-imperfect-match-default --args "AGCTAGCTAGCT -S 90 SRR000001" --rc 1 )
	add_test ( NAME SraSearch-6.1-imperfect-match-fgrep-standard  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
		    	${SRA_SEARCH} 6.1-imperfect-match-fgrep-standard --args "AGCTAGCTAGCT --score 90 -a FgrepStandard SRR000001" --rc 1  )
	add_test ( NAME SraSearch-6.2-imperfect-match-fgrep-boyer-moore  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
				${SRA_SEARCH} 6.2-imperfect-match-fgrep-boyer-moore --args "AGCTAGCTAGCT --score 90 --algorithm FgrepBoyerMoore SRR000001" --rc 1  )
	add_test ( NAME SraSearch-6.3-imperfect-match-fgrep-aho  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
	         	${SRA_SEARCH} 6.3-imperfect-match-fgrep-aho --args "AGCTAGCTAGCT --score 90 --algorithm FgrepAho SRR000001" --rc 1  )
	add_test ( NAME SraSearch-6.8-imperfect-match-nucstrstr  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
	         	${SRA_SEARCH} 6.8-imperfect-match-nucstrstr --args "AGCTAGCTAGCT --score 90 --algorithm NucStrstr SRR000001" --rc 1 )

	# multi-threading, bad arguments
	add_test ( NAME SraSearch-7.0.1-threads-no-arg WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
			 	${SRA_SEARCH} 7.0.1-threads-no-arg --rc 1 --args "AGCTAGCTAGCT SRR600094 SRR600095 -T" )
	add_test ( NAME SraSearch-7.0.2-threads-bad-arg WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
			 	${SRA_SEARCH} 7.0.2-threads-bad-arg --rc 1 --args "AGCTAGCTAGCT SRR600094 SRR600095 --threads -1" )

    # Slow tests (have to have "Slow" in the name, run without: "ctest -E Slow", run only: "ctest -R Slow"
	# single run, single threaded
	add_test ( NAME SraSearchSlow-3.0-single-run-default WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
				${SRA_SEARCH} 3.0-single-run-default            --args "AGCTAGCTAGCT SRR000001" )
	add_test ( NAME SraSearchSlow-3.1-single-run-fgrep-standard WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
				${SRA_SEARCH} 3.1-single-run-fgrep-standard     --args "AGCTAGCTAGCT -a FgrepStandard SRR000001" )
	add_test ( NAME SraSearchSlow-3.2-single-run-fgrep-boyer-moore WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
				${SRA_SEARCH} 3.2-single-run-fgrep-boyer-moore  --args "AGCTAGCTAGCT --algorithm FgrepBoyerMoore    SRR000001" )
	add_test ( NAME SraSearchSlow-3.3-single-run-fgrep-aho WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
				${SRA_SEARCH} 3.3-single-run-fgrep-aho         --args "AGCTAGCTAGCT --algorithm FgrepAho           SRR000001" )
	add_test ( NAME SraSearchSlow-3.4-single-run-agrep-dp WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
				${SRA_SEARCH} 3.4-single-run-agrep-dp           --args "AGCTAGCTAGCT --algorithm AgrepDP            SRR000001" )
	add_test ( NAME SraSearchSlow-3.5-single-run-agrep-wumanber WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
				${SRA_SEARCH} 3.5-single-run-agrep-wumanber     --args "AGCTAGCTAGCT --algorithm AgrepWuManber      SRR000001" )
	add_test ( NAME SraSearchSlow-3.6-single-run-agrep-myers WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
				${SRA_SEARCH} 3.6-single-run-agrep-myers        --args "AGCTAGCTAGCT --algorithm AgrepMyers         SRR000001" )
	add_test ( NAME SraSearchSlow-3.7-single-run-agrep-myers-unltd WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
				${SRA_SEARCH} 3.7-single-run-agrep-myers-unltd  --args "AGCTAGCTAGCT --algorithm AgrepMyersUnltd    SRR000001" )
	add_test ( NAME SraSearchSlow-3.8-single-run-nucstrstr WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
				${SRA_SEARCH} 3.8-single-run-nucstrstr          --args "AGCTAGCTAGCT --algorithm NucStrstr          SRR000001" )
	add_test ( NAME SraSearchSlow-3.9-single-run-smith-waterman WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
				${SRA_SEARCH} 3.9-single-run-smith-waterman     --args "AGCTAGCTAGCT --algorithm SmithWaterman      SRR000001" )

	# multiple runs, single threaded
	add_test ( NAME SraSearchSlow-4.0-multiple-runs WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
				${SRA_SEARCH} 4.0-multiple-runs --args "TACTGGGCGTAAAGCGTGCGCAGGCGGT --algorithm NucStrstr SRR600096 SRR600097" )

	# query expressions
	add_test ( NAME SraSearchSlow-5.0-run-nucstrstr-with-expression WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
				${SRA_SEARCH} 5.0-run-nucstrstr-with-expression --args "AAAAAAACCCCCCCAAAAAAACCCCCCC\\|\\|AGCTAGCTAGCT --algorithm NucStrstr --expression  SRR000001" )

	# imperfect match
	add_test ( NAME SraSearchSlow-6.4-imperfect-match-agrep-dp WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
				${SRA_SEARCH} 6.4-imperfect-match-agrep-dp           --args "AGCTAGCTAGCT --score 90 --algorithm AgrepDP            SRR000001" )
	add_test ( NAME SraSearchSlow-6.5-imperfect-match-agrep-wumanber WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
				${SRA_SEARCH} 6.5-imperfect-match-agrep-wumanber     --args "AGCTAGCTAGCT --score 90 --algorithm AgrepWuManber      SRR000001" )
	add_test ( NAME SraSearchSlow-6.6-imperfect-match-agrep-myers WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
				${SRA_SEARCH} 6.6-imperfect-match-agrep-myers        --args "AGCTAGCTAGCT --score 90 --algorithm AgrepMyers         SRR000001" )
	add_test ( NAME SraSearchSlow-6.7-imperfect-match-agrep-myers-unltd WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
				${SRA_SEARCH} 6.7-imperfect-match-agrep-myers-unltd  --args "AGCTAGCTAGCT --score 90 --algorithm AgrepMyersUnltd    SRR000001" )
	add_test ( NAME SraSearchSlow-6.9-imperfect-match-smith-waterman WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
				${SRA_SEARCH} 6.9-imperfect-match-smith-waterman     --args "AGCTAGCTAGCT --score 90 --algorithm SmithWaterman      SRR000001" )

	# multi-threading
	# multiple runs, single-threaded
	add_test ( NAME SraSearchSlow-7.0-nothreads WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
				${SRA_SEARCH} 7.0-nothreads --args "ACGTAGGGTCC --nothreads SRR600094 SRR600095 SRR600096 SRR600099" )
	# multiple runs, thread per accession
	add_test ( NAME SraSearchSlow-7.1-threads-unsorted WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
				${SRA_SEARCH} 7.1-threads-unsorted --sort --args "ACGTAGGGTCC --threadperacc --threads 4 SRR600094 SRR600095 SRR600096 SRR600099" )
	# single run, thread per blob
	add_test ( NAME SraSearchSlow-7.2-threads-single WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
				${SRA_SEARCH} 7.2-threads-single --sort --args "ACGTAGGGTCC --threads 2 SRR600094" )
	# multiple runs, thread per blob
	add_test ( NAME SraSearchSlow-7.3-threads-multiple WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
				${SRA_SEARCH} 7.3-threads-multiple --sort --args "ACGTAGGGTCC --threads 2 SRR600094 SRR600095 SRR600096 SRR600099" )

	# nothing found
	add_test ( NAME SraSearchSlow-8.0-nothing-found WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND runtestcase.sh
				${SRA_SEARCH} 8.0-nothing-found --rc 1 --args "TTTTTTTTAAAAAAACCCC SRR000001" )

endif ()
