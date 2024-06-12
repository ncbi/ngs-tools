if ( ( ${CMAKE_SYSTEM_NAME} MATCHES "Linux" )
        OR
     ( WIN32 ) AND NOT ( ${CMAKE_GENERATOR} MATCHES "Visual Studio 10 2010*" )
   )
   add_subdirectory(tax)
endif()

