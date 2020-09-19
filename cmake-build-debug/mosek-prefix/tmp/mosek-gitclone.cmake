
if(NOT "/home/griv/CLionProjects/MEB_TEST/cmake-build-debug/mosek-prefix/src/mosek-stamp/mosek-gitinfo.txt" IS_NEWER_THAN "/home/griv/CLionProjects/MEB_TEST/cmake-build-debug/mosek-prefix/src/mosek-stamp/mosek-gitclone-lastrun.txt")
  message(STATUS "Avoiding repeated git clone, stamp file is up to date: '/home/griv/CLionProjects/MEB_TEST/cmake-build-debug/mosek-prefix/src/mosek-stamp/mosek-gitclone-lastrun.txt'")
  return()
endif()

execute_process(
  COMMAND ${CMAKE_COMMAND} -E rm -rf "/home/griv/CLionProjects/MEB_TEST/cmake-build-debug/mosek-src"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to remove directory: '/home/griv/CLionProjects/MEB_TEST/cmake-build-debug/mosek-src'")
endif()

# try the clone 3 times in case there is an odd git clone issue
set(error_code 1)
set(number_of_tries 0)
while(error_code AND number_of_tries LESS 3)
  execute_process(
    COMMAND "/usr/bin/git"  clone --no-checkout "https://github.com/RobotLocomotion/mosek.git" "mosek-src"
    WORKING_DIRECTORY "/home/griv/CLionProjects/MEB_TEST/cmake-build-debug"
    RESULT_VARIABLE error_code
    )
  math(EXPR number_of_tries "${number_of_tries} + 1")
endwhile()
if(number_of_tries GREATER 1)
  message(STATUS "Had to git clone more than once:
          ${number_of_tries} times.")
endif()
if(error_code)
  message(FATAL_ERROR "Failed to clone repository: 'https://github.com/RobotLocomotion/mosek.git'")
endif()

execute_process(
  COMMAND "/usr/bin/git"  checkout e4bcf2234391c2a16adbce3063090e77e3b1f027 --
  WORKING_DIRECTORY "/home/griv/CLionProjects/MEB_TEST/cmake-build-debug/mosek-src"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to checkout tag: 'e4bcf2234391c2a16adbce3063090e77e3b1f027'")
endif()

set(init_submodules TRUE)
if(init_submodules)
  execute_process(
    COMMAND "/usr/bin/git"  submodule update --recursive --init 
    WORKING_DIRECTORY "/home/griv/CLionProjects/MEB_TEST/cmake-build-debug/mosek-src"
    RESULT_VARIABLE error_code
    )
endif()
if(error_code)
  message(FATAL_ERROR "Failed to update submodules in: '/home/griv/CLionProjects/MEB_TEST/cmake-build-debug/mosek-src'")
endif()

# Complete success, update the script-last-run stamp file:
#
execute_process(
  COMMAND ${CMAKE_COMMAND} -E copy
    "/home/griv/CLionProjects/MEB_TEST/cmake-build-debug/mosek-prefix/src/mosek-stamp/mosek-gitinfo.txt"
    "/home/griv/CLionProjects/MEB_TEST/cmake-build-debug/mosek-prefix/src/mosek-stamp/mosek-gitclone-lastrun.txt"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to copy script-last-run stamp file: '/home/griv/CLionProjects/MEB_TEST/cmake-build-debug/mosek-prefix/src/mosek-stamp/mosek-gitclone-lastrun.txt'")
endif()

