# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /snap/clion/126/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /snap/clion/126/bin/cmake/linux/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/griv/CLionProjects/MEB_TEST/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/griv/CLionProjects/MEB_TEST/cmake-build-debug/mosek_eigen_blaze-build

# Include any dependencies generated for this target.
include cxx/CMakeFiles/test_meb.dir/depend.make

# Include the progress variables for this target.
include cxx/CMakeFiles/test_meb.dir/progress.make

# Include the compile flags for this target's objects.
include cxx/CMakeFiles/test_meb.dir/flags.make

cxx/CMakeFiles/test_meb.dir/main.cpp.o: cxx/CMakeFiles/test_meb.dir/flags.make
cxx/CMakeFiles/test_meb.dir/main.cpp.o: /home/griv/CLionProjects/MEB_TEST/src/cxx/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/griv/CLionProjects/MEB_TEST/cmake-build-debug/mosek_eigen_blaze-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object cxx/CMakeFiles/test_meb.dir/main.cpp.o"
	cd /home/griv/CLionProjects/MEB_TEST/cmake-build-debug/mosek_eigen_blaze-build/cxx && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/test_meb.dir/main.cpp.o -c /home/griv/CLionProjects/MEB_TEST/src/cxx/main.cpp

cxx/CMakeFiles/test_meb.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_meb.dir/main.cpp.i"
	cd /home/griv/CLionProjects/MEB_TEST/cmake-build-debug/mosek_eigen_blaze-build/cxx && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/griv/CLionProjects/MEB_TEST/src/cxx/main.cpp > CMakeFiles/test_meb.dir/main.cpp.i

cxx/CMakeFiles/test_meb.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_meb.dir/main.cpp.s"
	cd /home/griv/CLionProjects/MEB_TEST/cmake-build-debug/mosek_eigen_blaze-build/cxx && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/griv/CLionProjects/MEB_TEST/src/cxx/main.cpp -o CMakeFiles/test_meb.dir/main.cpp.s

# Object files for target test_meb
test_meb_OBJECTS = \
"CMakeFiles/test_meb.dir/main.cpp.o"

# External object files for target test_meb
test_meb_EXTERNAL_OBJECTS =

cxx/test_meb: cxx/CMakeFiles/test_meb.dir/main.cpp.o
cxx/test_meb: cxx/CMakeFiles/test_meb.dir/build.make
cxx/test_meb: cxx/libmpc_mosek_eigen.so
cxx/test_meb: /home/griv/CLionProjects/MEB_TEST/cmake-build-debug/install/lib/libmosek64.so.7.1
cxx/test_meb: /home/griv/CLionProjects/MEB_TEST/cmake-build-debug/install/lib/libiomp5.so
cxx/test_meb: /usr/lib/x86_64-linux-gnu/libm.so
cxx/test_meb: /usr/lib/liblapack.so
cxx/test_meb: /usr/lib/libblas.so
cxx/test_meb: cxx/CMakeFiles/test_meb.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/griv/CLionProjects/MEB_TEST/cmake-build-debug/mosek_eigen_blaze-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable test_meb"
	cd /home/griv/CLionProjects/MEB_TEST/cmake-build-debug/mosek_eigen_blaze-build/cxx && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_meb.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
cxx/CMakeFiles/test_meb.dir/build: cxx/test_meb

.PHONY : cxx/CMakeFiles/test_meb.dir/build

cxx/CMakeFiles/test_meb.dir/clean:
	cd /home/griv/CLionProjects/MEB_TEST/cmake-build-debug/mosek_eigen_blaze-build/cxx && $(CMAKE_COMMAND) -P CMakeFiles/test_meb.dir/cmake_clean.cmake
.PHONY : cxx/CMakeFiles/test_meb.dir/clean

cxx/CMakeFiles/test_meb.dir/depend:
	cd /home/griv/CLionProjects/MEB_TEST/cmake-build-debug/mosek_eigen_blaze-build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/griv/CLionProjects/MEB_TEST/src /home/griv/CLionProjects/MEB_TEST/src/cxx /home/griv/CLionProjects/MEB_TEST/cmake-build-debug/mosek_eigen_blaze-build /home/griv/CLionProjects/MEB_TEST/cmake-build-debug/mosek_eigen_blaze-build/cxx /home/griv/CLionProjects/MEB_TEST/cmake-build-debug/mosek_eigen_blaze-build/cxx/CMakeFiles/test_meb.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : cxx/CMakeFiles/test_meb.dir/depend

