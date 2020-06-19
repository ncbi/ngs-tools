add_executable ( trace-loader
    tl_exception.cpp
    tl_owp.cpp
    tl_ows.cpp
    tl_report.cpp
    tl_trace_dp.cpp
    tl_trace_fwa.cpp
    tl_traceconfig.cpp
    tl_tracedata.cpp
    tl_tracedata_abi.cpp
    tl_tracedata_reader.cpp
    tl_tracedata_scf.cpp
    tl_tracedata_sff.cpp
    tl_tracedata_ztr.cpp
    tl_tracefields.cpp
    tl_tracefields_init.cpp
    tl_traceinfo.cpp
    tl_util.cpp
    tl_validator.cpp
    tl_vrules.cpp
    trace_loader.cpp
)

include_directories ( ${CMAKE_SOURCE_DIR}/libs/general-writer/inc )

add_dependencies ( trace-loader general-writer )

target_link_libraries ( trace-loader
    general-writer
    ncbi-wvdb-static
    kapp
    tk-version
	${SYS_LIBRARIES}
)

links_and_install ( trace-loader )

add_subdirectory ( test )
