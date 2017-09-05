#set(CMAKE_VERBOSE_MAKEFILE on)

add_executable        ( compute-coverage compute-coverage.cpp )
add_dependencies      ( compute-coverage ncbi-vdb )
target_link_libraries ( compute-coverage
	${SYS_LIBRARIES}
)
ncbi_copy_exec_to_old_outdir ( compute-coverage )
install              ( TARGETS compute-coverage RUNTIME DESTINATION bin )

add_executable        ( compute-coverage-ngs compute-coverage-ngs.cpp )
add_dependencies      ( compute-coverage-ngs ncbi-vdb )
target_link_libraries ( compute-coverage-ngs
	${SYS_LIBRARIES}
)
ncbi_copy_exec_to_old_outdir ( compute-coverage-ngs )
install              ( TARGETS compute-coverage-ngs RUNTIME DESTINATION bin )

