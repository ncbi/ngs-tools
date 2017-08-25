add_test ( PileupStats-Version
    sh -c "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/pileup-stats -V |
        grep --quiet `head -n 1 ${CMAKE_SOURCE_DIR}/libs/tk-version/toolkit.vers`"
    )

add_test ( PileupStats-IgnoreArgs-L        ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/pileup-stats -L 1 )
add_test ( PileupStats-IgnoreArgs-z        ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/pileup-stats -z ignore.txt )
add_test ( PileupStats-IgnoreArgs-loglovel ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/pileup-stats --log-level 1 )
add_test ( PileupStats-IgnoreArgs-xmllog   ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/pileup-stats --xml-log ignore.txt )
add_test ( PileupStats-IgnoreArgs-xmllogfd ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/pileup-stats --xml-log-fd ignore.txt )
