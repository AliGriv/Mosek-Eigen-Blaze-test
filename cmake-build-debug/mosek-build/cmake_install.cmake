# Install script for directory: /home/griv/CLionProjects/MEB_TEST/cmake-build-debug/mosek-src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/griv/CLionProjects/MEB_TEST/cmake-build-debug/install")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mosek" TYPE FILE FILES "/home/griv/CLionProjects/MEB_TEST/cmake-build-debug/mosek-build/mosek/7/tools/platform/linux64x86/h/mosek.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/mosek" TYPE FILE RENAME "mosek-c-targets.cmake" FILES "/home/griv/CLionProjects/MEB_TEST/cmake-build-debug/mosek-build/mosek-c-targets-install.cmake")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig" TYPE FILE FILES "/home/griv/CLionProjects/MEB_TEST/cmake-build-debug/mosek-build/mosek.pc")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE FILE FILES
    "/home/griv/CLionProjects/MEB_TEST/cmake-build-debug/mosek-build/mosek/7/tools/platform/linux64x86/bin/libiomp5.so"
    "/home/griv/CLionProjects/MEB_TEST/cmake-build-debug/mosek-build/mosek/7/tools/platform/linux64x86/bin/libmosek64.so"
    "/home/griv/CLionProjects/MEB_TEST/cmake-build-debug/mosek-build/mosek/7/tools/platform/linux64x86/bin/libmosek64.so.7.1"
    "/home/griv/CLionProjects/MEB_TEST/cmake-build-debug/mosek-build/mosek/7/tools/platform/linux64x86/bin/libmosekglb64.so.7.1"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/doc/mosek" TYPE FILE FILES "/home/griv/CLionProjects/MEB_TEST/cmake-build-debug/mosek-build/mosek/7/license.pdf")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/mosek" TYPE FILE FILES
    "/home/griv/CLionProjects/MEB_TEST/cmake-build-debug/mosek-build/mosek-config.cmake"
    "/home/griv/CLionProjects/MEB_TEST/cmake-build-debug/mosek-build/mosek-config-version.cmake"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/griv/CLionProjects/MEB_TEST/cmake-build-debug/mosek-build/test/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/griv/CLionProjects/MEB_TEST/cmake-build-debug/mosek-build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
