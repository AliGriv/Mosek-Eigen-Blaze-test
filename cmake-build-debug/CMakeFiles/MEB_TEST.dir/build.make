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
CMAKE_SOURCE_DIR = /home/griv/CLionProjects/MEB_TEST

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/griv/CLionProjects/MEB_TEST/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/MEB_TEST.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/MEB_TEST.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/MEB_TEST.dir/flags.make

CMakeFiles/MEB_TEST.dir/main.cpp.o: CMakeFiles/MEB_TEST.dir/flags.make
CMakeFiles/MEB_TEST.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/griv/CLionProjects/MEB_TEST/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/MEB_TEST.dir/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/MEB_TEST.dir/main.cpp.o -c /home/griv/CLionProjects/MEB_TEST/main.cpp

CMakeFiles/MEB_TEST.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MEB_TEST.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/griv/CLionProjects/MEB_TEST/main.cpp > CMakeFiles/MEB_TEST.dir/main.cpp.i

CMakeFiles/MEB_TEST.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MEB_TEST.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/griv/CLionProjects/MEB_TEST/main.cpp -o CMakeFiles/MEB_TEST.dir/main.cpp.s

# Object files for target MEB_TEST
MEB_TEST_OBJECTS = \
"CMakeFiles/MEB_TEST.dir/main.cpp.o"

# External object files for target MEB_TEST
MEB_TEST_EXTERNAL_OBJECTS =

MEB_TEST: CMakeFiles/MEB_TEST.dir/main.cpp.o
MEB_TEST: CMakeFiles/MEB_TEST.dir/build.make
MEB_TEST: CMakeFiles/MEB_TEST.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/griv/CLionProjects/MEB_TEST/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable MEB_TEST"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/MEB_TEST.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/MEB_TEST.dir/build: MEB_TEST

.PHONY : CMakeFiles/MEB_TEST.dir/build

CMakeFiles/MEB_TEST.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/MEB_TEST.dir/cmake_clean.cmake
.PHONY : CMakeFiles/MEB_TEST.dir/clean

CMakeFiles/MEB_TEST.dir/depend:
	cd /home/griv/CLionProjects/MEB_TEST/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/griv/CLionProjects/MEB_TEST /home/griv/CLionProjects/MEB_TEST /home/griv/CLionProjects/MEB_TEST/cmake-build-debug /home/griv/CLionProjects/MEB_TEST/cmake-build-debug /home/griv/CLionProjects/MEB_TEST/cmake-build-debug/CMakeFiles/MEB_TEST.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/MEB_TEST.dir/depend

