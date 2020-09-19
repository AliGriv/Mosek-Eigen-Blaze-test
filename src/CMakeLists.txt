cmake_minimum_required(VERSION 3.17)
project(mosek_eigen_blaze)

if(APPLE)
    # Compiler ID for Apple Clang is now AppleClang.
    if(POLICY CMP0025)
        cmake_policy(SET CMP0025 NEW)
    endif()

    # MACOSX_RPATH is enabled by default.
    if(POLICY CMP0042)
        cmake_policy(SET CMP0042 NEW)
    endif()
endif()

include(CMakePackageConfigHelpers)
include(CTest)
include(GNUInstallDirs)

list(INSERT CMAKE_MODULE_PATH 0 "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")
find_package(Mosek CONFIG REQUIRED)

find_package(Eigen3 CONFIG)

if(NOT Eigen3_FOUND)
    find_package(Eigen3 MODULE REQUIRED)
endif()

# Set RPATH for installed binaries.
set(CMAKE_INSTALL_RPATH
        "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"
        "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)

list(REMOVE_DUPLICATES CMAKE_INSTALL_RPATH)

# If possible, use the functionality in CMake 3.1 and above (CMake 3.9 and
# above for the PGI and XL compilers) to set the compiler flag for setting the
# C++ dialect to C++0x or C++11.
if(CMAKE_VERSION VERSION_LESS 3.1
        OR (CMAKE_VERSION VERSION_LESS 3.9 AND CMAKE_CXX_COMPILER_ID MATCHES "^(PGI|XL)$"))
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_FLAGS_CXX14}")
else()
    set(CMAKE_CXX_STANDARD 14)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_FLAGS_WARNING}")
set(MEB_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
set(MEB_INCLUDE_DIR ${CMAKE_INSTALL_INCLUDEDIR})
set(MEB_LIBRARY_DIR ${CMAKE_INSTALL_LIBDIR})
add_subdirectory(cxx)


