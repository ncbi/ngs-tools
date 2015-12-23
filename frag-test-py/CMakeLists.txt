cmake_minimum_required (VERSION 2.8)

find_package(PythonLibs)

project ( ngs-fragtest-py )

add_dependencies ( ngs-fragtest-py ngs ncbi-vdb ) 

install ( FILES FragTest.py DESTINATION share/ngs-tools/examples/python/ )

