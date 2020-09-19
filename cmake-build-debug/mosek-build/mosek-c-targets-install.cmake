
####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was mosek-c-targets.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

####################################################################################

if(NOT TARGET Threads::Threads)
  message(FATAL_ERROR
    "Could NOT find the target Threads::Threads referenced by target mosek"
  )
endif()

if(NOT TARGET mosek)
  add_library(mosek SHARED IMPORTED)
  set_target_properties(mosek PROPERTIES
    INTERFACE_LINK_LIBRARIES
      "${PACKAGE_PREFIX_DIR}/lib/libiomp5.so;/usr/lib/x86_64-linux-gnu/libm.so;Threads::Threads"
    IMPORTED_LOCATION "${PACKAGE_PREFIX_DIR}/lib/libmosek64.so.7.1"
    INTERFACE_INCLUDE_DIRECTORIES "${PACKAGE_PREFIX_DIR}/include/mosek"
  )
endif()
