
####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was mosek_eigen_blaze-config.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

####################################################################################

include(CMakeFindDependencyMacro)


find_dependency(mosek)

set(MEB_FOUND TRUE)
set_and_check(MEB_INCLUDE_DIRS "${PACKAGE_PREFIX_DIR}/include")
set_and_check(MEB_LIBRARY_DIRS "${PACKAGE_PREFIX_DIR}/lib")
set(MEB_LIBRARIES)

include("${CMAKE_CURRENT_LIST_DIR}/mosek_eigen_blaze-targets.cmake")
