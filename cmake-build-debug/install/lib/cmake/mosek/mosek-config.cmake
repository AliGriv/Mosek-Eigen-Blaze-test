cmake_minimum_required(VERSION 3.1)


####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was mosek-config.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

####################################################################################

include(CMakeFindDependencyMacro)

set(MOSEK_VERSION_MAJOR "7")
set(MOSEK_VERSION_MINOR "1")
set(MOSEK_VERSION_MINOR "0")
set(MOSEK_VERSION_TWEAK "45")
set(MOSEK_VERSION "7.1.0.45")

find_dependency(Threads)

set_and_check(MOSEK_INCLUDE_DIRS "${PACKAGE_PREFIX_DIR}/include/mosek")
set_and_check(MOSEK_LIBRARY_DIRS "${PACKAGE_PREFIX_DIR}/lib")
set(MOSEK_LIBRARIES mosek "/usr/lib/x86_64-linux-gnu/libm.so" Threads::Threads)

include("${CMAKE_CURRENT_LIST_DIR}/mosek-c-targets.cmake")

