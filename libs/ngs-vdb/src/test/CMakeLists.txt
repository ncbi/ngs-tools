include_directories ( ${CMAKE_SOURCE_DIR} )
include_directories ("${VDB_INCDIR}/../libs/ngs")

set ( TEST_LIBS
    ngs-vdb
    ${CMAKE_STATIC_LIBRARY_PREFIX}ktst${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${CMAKE_STATIC_LIBRARY_PREFIX}kapp${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${SYS_LIBRARIES}
)

if (WIN32)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /SUBSYSTEM:WINDOWS /ENTRY:wmainCRTStartup")
endif ()

#-------------------------------------------------------------------------------
# tests for the underlying C APIs

add_definitions(-D__mod__="NGS_TEST")
add_definitions(-D__file__="__FILE__")
add_definitions(-D__fext__="cpp")

add_executable ( test-ngs ngstest.cpp )
target_link_libraries ( test-ngs ${TEST_LIBS})
add_test ( NAME Test_Ngs COMMAND test-ngs WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable ( test-ngs_sra ngstest_sra.cpp )
target_link_libraries ( test-ngs_sra ${TEST_LIBS})
add_test ( NAME Test_Ngs_Sra COMMAND test-ngs_sra WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable ( test-ngs_sradb ngstest_sradb.cpp )
target_link_libraries ( test-ngs_sradb ${TEST_LIBS})
add_test ( NAME Test_Ngs_SraDb COMMAND test-ngs_sradb WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable ( test-ngs_csra1 ngstest_csra1.cpp )
target_link_libraries ( test-ngs_csra1 ${TEST_LIBS})
add_test ( NAME SlowTest_Ngs_Csra1 COMMAND test-ngs_csra1 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable ( test-ngs_csra1_readcollection ngstest_csra1_readcollection.cpp )
target_link_libraries ( test-ngs_csra1_readcollection ${TEST_LIBS})
add_test ( NAME Test_Ngs_Csra1_ReadCollection COMMAND test-ngs_csra1_readcollection WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable ( test-ngs_csra1_refwin ngstest_csra1_refwin.cpp )
target_link_libraries ( test-ngs_csra1_refwin ${TEST_LIBS})
add_test ( NAME SlowTest_Ngs_Csra1_RefWin COMMAND test-ngs_csra1_refwin WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable ( test-ngs_csra1_pileup ngstest_csra1_pileup.cpp )
target_link_libraries ( test-ngs_csra1_pileup ${TEST_LIBS})
add_test ( NAME SlowTest_Ngs_Csra1_Pileup COMMAND test-ngs_csra1_pileup WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable ( test-ngs_reference ngstest_reference.cpp )
target_link_libraries ( test-ngs_reference ${TEST_LIBS})
add_test ( NAME SlowTest_Ngs_Reference COMMAND test-ngs_reference WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable ( test-fragmentblob ngstest_fragmentblob.cpp )
target_link_libraries ( test-fragmentblob ${TEST_LIBS})
add_test ( NAME SlowTest_Ngs_FragmentBlob COMMAND test-fragmentblob WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable ( test-referenceblob ngstest_referenceblob.cpp )
target_link_libraries ( test-referenceblob ${TEST_LIBS})
add_test ( NAME Test_Ngs_ReferenceBlob COMMAND test-referenceblob WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable ( test-byteblob ngstest_byteblob.cpp )
target_link_libraries ( test-byteblob ${TEST_LIBS})
add_test ( NAME Test_Ngs_ByteBlob COMMAND test-byteblob WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

#-------------------------------------------------------------------------------
# tests for the NGS C++ APIs

add_executable ( test-ngs_sra_c++ SRAtest.cpp )
target_link_libraries ( test-ngs_sra_c++ ${TEST_LIBS})
add_test ( NAME Test_Ngs_CPP_Sra COMMAND test-ngs_sra_c++ WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable ( test-ngs_sradb_c++ SRADBtest.cpp )
target_link_libraries ( test-ngs_sradb_c++ ${TEST_LIBS})
add_test ( NAME Test_Ngs_CPP_SraDb COMMAND test-ngs_sradb_c++ WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable ( test-ngs_csra1_c++ CSRA1test.cpp )
target_link_libraries ( test-ngs_csra1_c++ ${TEST_LIBS})
add_test ( NAME SlowTest_Ngs_CPP_CSRA1 COMMAND test-ngs_csra1_c++ WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable ( test-ngs_csra1pileup++ CSRA1PileupTest.cpp )
target_link_libraries ( test-ngs_csra1pileup++ ${TEST_LIBS})
add_test ( NAME SlowTest_Ngs_CPP_CSRA1Pileup COMMAND test-ngs_csra1pileup++ WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

#-------------------------------------------------------------------------------
# tests for the NGS Java APIs

set( JUNIT ${CMAKE_SOURCE_DIR}/shared/java/junit-4.11.jar)
set( CLASSPATH ${JUNIT}:${NGSJAR})
set( LIBRARY_PATH -Djava.library.path=${NGS_LIBDIR}:${VDB_LIBDIR})

# use standard java System.loadLibrary
set ( VARS ${VARS} -Dvdb.System.loadLibrary=1 )

# remove all found JNI DLL-s
#set ( VARS ${VARS} -Dvdb.deleteLibraries=1 )

# custom load library search
#set ( VARS ${VARS} -Dvdb.loadLibraryLocations=W )

set ( VARS ${VARS} -Dvdb.log=FINEST )

#set(CMAKE_JAVA_COMPILE_FLAGS -verbose)

include (UseJava)
add_jar( ngs_test_lib_load
    SOURCES ngs_test_lib_load.java
    INCLUDE_JARS ${JUNIT} ${NGSJAR}
    OUTPUT_DIR ${PROJECT_BINARY_DIR}
)

add_test ( NAME Test_Ngs_Java_lib_load
    COMMAND java -cp ${PROJECT_BINARY_DIR}/ngs_test_lib_load.jar:${CLASSPATH} ${VARS} ${LIBRARY_PATH} org.junit.runner.JUnitCore gov.nih.nlm.ncbi.ngs.ngs_test_lib_load
)

add_jar( ngs_test_CSRA1
    SOURCES ngs_test_CSRA1.java
    INCLUDE_JARS ${JUNIT} ${NGSJAR}
    OUTPUT_DIR ${PROJECT_BINARY_DIR}
)

add_test ( NAME SlowTest_Ngs_Java_CSRA1
    COMMAND java -cp ${PROJECT_BINARY_DIR}/ngs_test_CSRA1.jar:${CLASSPATH} ${VARS} ${LIBRARY_PATH} org.junit.runner.JUnitCore ngs_test_CSRA1
)

#-------------------------------------------------------------------------------
# tests for this library, ngs-vdb

add_executable ( test-ngs-vdb test-ngs-vdb.cpp )
add_dependencies ( test-ngs-vdb ngs-vdb )

target_link_libraries ( test-ngs-vdb
    tk-version
    ${TEST_LIBS}
)

add_test ( NAME SlowTest_NgsVdbWhiteBox COMMAND test-ngs-vdb )

