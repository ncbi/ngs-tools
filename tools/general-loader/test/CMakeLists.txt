include_directories (".." ${CMAKE_SOURCE_DIR}/libs )

add_executable ( test-general-loader
    test-general-loader
    testsource
)

add_dependencies ( test-general-loader general-writer )

target_link_libraries ( test-general-loader
    general-writer
    kapp
    ncbi-wvdb-static
    ktst
    load
    dl
    pthread
)

add_executable ( test-gw-dumper test-gw-dumper )
target_link_libraries ( test-gw-dumper
    general-writer
)

#white box tests
file ( MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/db )
file ( MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/input )
add_test ( NAME GeneralLoaderWhiteBox
           COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test-general-loader
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
         )
add_test ( GWDumper ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test-gw-dumper )

# command line tests
add_executable ( makeinputs
    makeinputs
    testsource
)
target_link_libraries ( makeinputs
    general-writer
    ncbi-vdb-static
)

if ( NOT WIN32)

    # tests depending on vdb-dump
    find_program(VDB_DUMP vdb-dump HINTS ${SRATOOLS_BINDIR})
    if (NOT VDB_DUMP)
        message(WARNING "vdb-dump is not found, skipping tests depending on it")
    else ()
        message(STATUS "vdb-dump found in ${VDB_DUMP}")

        set ( CMD runtestcase.sh ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} ${VDB_DUMP} ${CMAKE_CURRENT_SOURCE_DIR} )

        add_test ( NAME makeinputs
                COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/makeinputs
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )

        add_test ( NAME GeneralLoader-1.0-Basic-Unpacked
                COMMAND ${CMD} 1 0 -I${VDB_INCDIR}
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )
        add_test ( NAME GeneralLoader-1.1-Basic-Packed
                COMMAND ${CMD} 1packed 0 -I${VDB_INCDIR}
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )
        add_test ( NAME GeneralLoader-1.2-TargetDbOverride
                COMMAND ${CMD} 1override 0 "-I${VDB_INCDIR} --target actual/1override/db"
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )
        add_test ( NAME GeneralLoader-1.3-TargetDbOverrideShorthand
                COMMAND ${CMD} 1override 0 "-I${VDB_INCDIR} -T actual/1override/db"
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )

        add_test ( NAME GeneralLoader-2.0-ErrorMessageEvent-Unpacked
                COMMAND ${CMD} 2 3 "-I${VDB_INCDIR} -L=err"
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )
        add_test ( NAME GeneralLoader-2.1-ErrorMessageEvent-Packed
                COMMAND ${CMD} 2packed 3 "-I${VDB_INCDIR} -L=err"
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )

        add_test ( NAME GeneralLoader-3.0-EmptyDefaultValues-Unpacked
                COMMAND ${CMD} 3 0 -I${VDB_INCDIR}
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )
        add_test ( NAME GeneralLoader-3.1-EmptyDefaultValues-Packed
                COMMAND ${CMD} 3 0 -I${VDB_INCDIR}
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )

        add_test ( NAME GeneralLoader-4.0-MoveAheadEvent-Unpacked
                COMMAND ${CMD} 4 0 -I${VDB_INCDIR}
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )
        add_test ( NAME GeneralLoader-4.1-MoveAheadEvent-Packed
                COMMAND ${CMD} 4packed 0 -I${VDB_INCDIR}
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )

        add_test ( NAME GeneralLoader-5.0-IntegerConversion #packed only
                COMMAND ${CMD} 5packed 0 -I${VDB_INCDIR}
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )

        set_tests_properties(
            GeneralLoader-1.0-Basic-Unpacked
            GeneralLoader-1.1-Basic-Packed
            GeneralLoader-1.2-TargetDbOverride
            GeneralLoader-1.3-TargetDbOverrideShorthand
            GeneralLoader-2.0-ErrorMessageEvent-Unpacked
            GeneralLoader-2.1-ErrorMessageEvent-Packed
            GeneralLoader-3.0-EmptyDefaultValues-Unpacked
            GeneralLoader-3.1-EmptyDefaultValues-Packed
            GeneralLoader-4.0-MoveAheadEvent-Unpacked
            GeneralLoader-4.1-MoveAheadEvent-Packed
            GeneralLoader-5.0-IntegerConversion

            PROPERTIES DEPENDS makeinputs
        )
    endif()

    # tests depending on kdbmeta
    find_program(KDBMETA kdbmeta HINTS ${SRATOOLS_BINDIR})
    if (NOT KDBMETA)
        message(WARNING "kdbmeta is not found, skipping tests depending on it")
    else ()
        message(STATUS "kdbmeta found in ${KDBMETA}")

        set ( CMD runtestcase.sh ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} ${KDBMETA} ${CMAKE_CURRENT_SOURCE_DIR} )

        add_test ( NAME GeneralLoader-6.0-SoftwareName-Unpacked
                COMMAND ${CMD} 6 0 "-L=info -I${VDB_INCDIR}" SOFTWARE/formatter
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )
        add_test ( NAME GeneralLoader-6.1-SoftwareName-Packed
                COMMAND ${CMD} 6packed 0 "-L=info -I${VDB_INCDIR}" SOFTWARE/formatter
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )

        set_tests_properties(
            GeneralLoader-6.0-SoftwareName-Unpacked
            GeneralLoader-6.1-SoftwareName-Packed

            PROPERTIES DEPENDS makeinputs
        )
    endif()

endif()

add_subdirectory(vschema)
