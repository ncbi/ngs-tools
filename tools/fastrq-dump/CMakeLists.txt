add_executable ( fastrq-dump fastrq-dump.cpp filter.cpp formatter.cpp operation.cpp settings.cpp )

target_link_libraries ( fastrq-dump
	${SYS_LIBRARIES}
)

ncbi_copy_exec_to_old_outdir ( fastrq-dump )
install ( TARGETS fastrq-dump RUNTIME DESTINATION bin )
