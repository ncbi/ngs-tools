add_executable ( dump-ref-fasta DumpReferenceFASTA.cpp)

target_link_libraries ( dump-ref-fasta
	${SYS_LIBRARIES}
)

links_and_install (dump-ref-fasta)
