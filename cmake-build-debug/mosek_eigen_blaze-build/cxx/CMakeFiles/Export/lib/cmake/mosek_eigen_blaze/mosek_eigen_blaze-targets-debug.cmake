#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "mosek_eigen_blaze" for configuration "Debug"
set_property(TARGET mosek_eigen_blaze APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(mosek_eigen_blaze PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libmosek_eigen_blaze.so"
  IMPORTED_SONAME_DEBUG "libmosek_eigen_blaze.so"
  )

list(APPEND _IMPORT_CHECK_TARGETS mosek_eigen_blaze )
list(APPEND _IMPORT_CHECK_FILES_FOR_mosek_eigen_blaze "${_IMPORT_PREFIX}/lib/libmosek_eigen_blaze.so" )

# Import target "mpc_mosek_eigen" for configuration "Debug"
set_property(TARGET mpc_mosek_eigen APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(mpc_mosek_eigen PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libmpc_mosek_eigen.so"
  IMPORTED_SONAME_DEBUG "libmpc_mosek_eigen.so"
  )

list(APPEND _IMPORT_CHECK_TARGETS mpc_mosek_eigen )
list(APPEND _IMPORT_CHECK_FILES_FOR_mpc_mosek_eigen "${_IMPORT_PREFIX}/lib/libmpc_mosek_eigen.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
