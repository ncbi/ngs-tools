cmake_minimum_required (VERSION 2.8)

ADD_CUSTOM_TARGET( ant_build ALL )
ADD_CUSTOM_COMMAND(
  TARGET  ant_build
  COMMAND ant -Dplatforms.JDK_1.7.home=$(JAVA_HOME) -buildfile ${CMAKE_CURRENT_SOURCE_DIR}/build.xml -Ddist.dir=${CMAKE_CURRENT_BINARY_DIR}
)
