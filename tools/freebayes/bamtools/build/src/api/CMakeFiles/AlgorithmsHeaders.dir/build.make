# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/gosuzombie/Desktop/pipeline/tools/freebayes/bamtools

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/gosuzombie/Desktop/pipeline/tools/freebayes/bamtools/build

# Utility rule file for AlgorithmsHeaders.

# Include the progress variables for this target.
include src/api/CMakeFiles/AlgorithmsHeaders.dir/progress.make

src/api/CMakeFiles/AlgorithmsHeaders:
	$(CMAKE_COMMAND) -E cmake_progress_report /home/gosuzombie/Desktop/pipeline/tools/freebayes/bamtools/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Exporting AlgorithmsHeaders"

AlgorithmsHeaders: src/api/CMakeFiles/AlgorithmsHeaders
AlgorithmsHeaders: src/api/CMakeFiles/AlgorithmsHeaders.dir/build.make
	cd /home/gosuzombie/Desktop/pipeline/tools/freebayes/bamtools/build/src/api && /usr/bin/cmake -E copy_if_different /home/gosuzombie/Desktop/pipeline/tools/freebayes/bamtools/src/api/algorithms/Sort.h /home/gosuzombie/Desktop/pipeline/tools/freebayes/bamtools/include/api/algorithms/Sort.h
.PHONY : AlgorithmsHeaders

# Rule to build all files generated by this target.
src/api/CMakeFiles/AlgorithmsHeaders.dir/build: AlgorithmsHeaders
.PHONY : src/api/CMakeFiles/AlgorithmsHeaders.dir/build

src/api/CMakeFiles/AlgorithmsHeaders.dir/clean:
	cd /home/gosuzombie/Desktop/pipeline/tools/freebayes/bamtools/build/src/api && $(CMAKE_COMMAND) -P CMakeFiles/AlgorithmsHeaders.dir/cmake_clean.cmake
.PHONY : src/api/CMakeFiles/AlgorithmsHeaders.dir/clean

src/api/CMakeFiles/AlgorithmsHeaders.dir/depend:
	cd /home/gosuzombie/Desktop/pipeline/tools/freebayes/bamtools/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/gosuzombie/Desktop/pipeline/tools/freebayes/bamtools /home/gosuzombie/Desktop/pipeline/tools/freebayes/bamtools/src/api /home/gosuzombie/Desktop/pipeline/tools/freebayes/bamtools/build /home/gosuzombie/Desktop/pipeline/tools/freebayes/bamtools/build/src/api /home/gosuzombie/Desktop/pipeline/tools/freebayes/bamtools/build/src/api/CMakeFiles/AlgorithmsHeaders.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/api/CMakeFiles/AlgorithmsHeaders.dir/depend

