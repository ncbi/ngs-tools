add_executable ( fastrq-dump fastrq-dump.cpp filter.cpp formatter.cpp operation.cpp settings.cpp )

target_link_libraries ( fastrq-dump
	${SYS_LIBRARIES}
)

install ( TARGETS fastrq-dump RUNTIME DESTINATION bin )
