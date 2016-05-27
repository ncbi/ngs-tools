cmake_minimum_required (VERSION 2.8)

find_package(Java)

ADD_CUSTOM_TARGET( ant_build ALL )
ADD_CUSTOM_COMMAND(
  TARGET  ant_build
  COMMAND ant   -Dplatform.java="${Java_JAVA_EXECUTABLE}"
                -Dplatform.javac="${Java_JAVAC_EXECUTABLE}"
                -Dplatform.javadoc="${Java_JAVADOC_EXECUTABLE}"
                -buildfile ${CMAKE_CURRENT_SOURCE_DIR}/build.xml 
                -Ddist.dir=${CMAKE_CURRENT_BINARY_DIR}
)

#-Dplatforms.JDK_1.7.home="$(JAVA_HOME)" 