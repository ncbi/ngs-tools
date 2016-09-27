add_executable ( dump-ref-fasta DumpReferenceFASTA.cpp )

add_dependencies ( dump-ref-fasta ncbi-vdb )

target_link_libraries ( dump-ref-fasta
	${SYS_LIBRARIES}
)

ncbi_copy_exec_to_old_outdir ( dump-ref-fasta )
install ( TARGETS dump-ref-fasta RUNTIME DESTINATION bin )

