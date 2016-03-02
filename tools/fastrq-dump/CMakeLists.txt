cmake_minimum_required (VERSION 2.8)

project ( fastrq-dump )

add_executable ( fastrq-dump fastrq-dump.cpp filter.cpp formatter.cpp operation.cpp settings.cpp )

add_dependencies ( fastrq-dump ngs ncbi-vdb ) 

target_link_libraries ( fastrq-dump 
	${SYS_LIBRARIES}
)

install ( TARGETS fastrq-dump RUNTIME DESTINATION bin )

