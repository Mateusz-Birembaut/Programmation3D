# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.31

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

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files\CMake\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\CMake\bin\cmake.exe" -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = C:\Users\mateu\Documents\Bureau\prog_3d\Prog3D\Cameras_TP_21_11_2024

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = C:\Users\mateu\Documents\Bureau\prog_3d\Prog3D\Cameras_TP_21_11_2024\TP

# Include any dependencies generated for this target.
include external/glfw-3.1.2/tests/CMakeFiles/msaa.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include external/glfw-3.1.2/tests/CMakeFiles/msaa.dir/compiler_depend.make

# Include the progress variables for this target.
include external/glfw-3.1.2/tests/CMakeFiles/msaa.dir/progress.make

# Include the compile flags for this target's objects.
include external/glfw-3.1.2/tests/CMakeFiles/msaa.dir/flags.make

external/glfw-3.1.2/tests/CMakeFiles/msaa.dir/codegen:
.PHONY : external/glfw-3.1.2/tests/CMakeFiles/msaa.dir/codegen

external/glfw-3.1.2/tests/CMakeFiles/msaa.dir/msaa.c.obj: external/glfw-3.1.2/tests/CMakeFiles/msaa.dir/flags.make
external/glfw-3.1.2/tests/CMakeFiles/msaa.dir/msaa.c.obj: external/glfw-3.1.2/tests/CMakeFiles/msaa.dir/includes_C.rsp
external/glfw-3.1.2/tests/CMakeFiles/msaa.dir/msaa.c.obj: C:/Users/mateu/Documents/Bureau/prog_3d/Prog3D/Cameras_TP_21_11_2024/external/glfw-3.1.2/tests/msaa.c
external/glfw-3.1.2/tests/CMakeFiles/msaa.dir/msaa.c.obj: external/glfw-3.1.2/tests/CMakeFiles/msaa.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=C:\Users\mateu\Documents\Bureau\prog_3d\Prog3D\Cameras_TP_21_11_2024\TP\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object external/glfw-3.1.2/tests/CMakeFiles/msaa.dir/msaa.c.obj"
	cd /d C:\Users\mateu\Documents\Bureau\prog_3d\Prog3D\Cameras_TP_21_11_2024\TP\external\glfw-3.1.2\tests && C:\msys64\ucrt64\bin\cc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT external/glfw-3.1.2/tests/CMakeFiles/msaa.dir/msaa.c.obj -MF CMakeFiles\msaa.dir\msaa.c.obj.d -o CMakeFiles\msaa.dir\msaa.c.obj -c C:\Users\mateu\Documents\Bureau\prog_3d\Prog3D\Cameras_TP_21_11_2024\external\glfw-3.1.2\tests\msaa.c

external/glfw-3.1.2/tests/CMakeFiles/msaa.dir/msaa.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/msaa.dir/msaa.c.i"
	cd /d C:\Users\mateu\Documents\Bureau\prog_3d\Prog3D\Cameras_TP_21_11_2024\TP\external\glfw-3.1.2\tests && C:\msys64\ucrt64\bin\cc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E C:\Users\mateu\Documents\Bureau\prog_3d\Prog3D\Cameras_TP_21_11_2024\external\glfw-3.1.2\tests\msaa.c > CMakeFiles\msaa.dir\msaa.c.i

external/glfw-3.1.2/tests/CMakeFiles/msaa.dir/msaa.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/msaa.dir/msaa.c.s"
	cd /d C:\Users\mateu\Documents\Bureau\prog_3d\Prog3D\Cameras_TP_21_11_2024\TP\external\glfw-3.1.2\tests && C:\msys64\ucrt64\bin\cc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S C:\Users\mateu\Documents\Bureau\prog_3d\Prog3D\Cameras_TP_21_11_2024\external\glfw-3.1.2\tests\msaa.c -o CMakeFiles\msaa.dir\msaa.c.s

external/glfw-3.1.2/tests/CMakeFiles/msaa.dir/__/deps/getopt.c.obj: external/glfw-3.1.2/tests/CMakeFiles/msaa.dir/flags.make
external/glfw-3.1.2/tests/CMakeFiles/msaa.dir/__/deps/getopt.c.obj: external/glfw-3.1.2/tests/CMakeFiles/msaa.dir/includes_C.rsp
external/glfw-3.1.2/tests/CMakeFiles/msaa.dir/__/deps/getopt.c.obj: C:/Users/mateu/Documents/Bureau/prog_3d/Prog3D/Cameras_TP_21_11_2024/external/glfw-3.1.2/deps/getopt.c
external/glfw-3.1.2/tests/CMakeFiles/msaa.dir/__/deps/getopt.c.obj: external/glfw-3.1.2/tests/CMakeFiles/msaa.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=C:\Users\mateu\Documents\Bureau\prog_3d\Prog3D\Cameras_TP_21_11_2024\TP\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object external/glfw-3.1.2/tests/CMakeFiles/msaa.dir/__/deps/getopt.c.obj"
	cd /d C:\Users\mateu\Documents\Bureau\prog_3d\Prog3D\Cameras_TP_21_11_2024\TP\external\glfw-3.1.2\tests && C:\msys64\ucrt64\bin\cc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT external/glfw-3.1.2/tests/CMakeFiles/msaa.dir/__/deps/getopt.c.obj -MF CMakeFiles\msaa.dir\__\deps\getopt.c.obj.d -o CMakeFiles\msaa.dir\__\deps\getopt.c.obj -c C:\Users\mateu\Documents\Bureau\prog_3d\Prog3D\Cameras_TP_21_11_2024\external\glfw-3.1.2\deps\getopt.c

external/glfw-3.1.2/tests/CMakeFiles/msaa.dir/__/deps/getopt.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/msaa.dir/__/deps/getopt.c.i"
	cd /d C:\Users\mateu\Documents\Bureau\prog_3d\Prog3D\Cameras_TP_21_11_2024\TP\external\glfw-3.1.2\tests && C:\msys64\ucrt64\bin\cc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E C:\Users\mateu\Documents\Bureau\prog_3d\Prog3D\Cameras_TP_21_11_2024\external\glfw-3.1.2\deps\getopt.c > CMakeFiles\msaa.dir\__\deps\getopt.c.i

external/glfw-3.1.2/tests/CMakeFiles/msaa.dir/__/deps/getopt.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/msaa.dir/__/deps/getopt.c.s"
	cd /d C:\Users\mateu\Documents\Bureau\prog_3d\Prog3D\Cameras_TP_21_11_2024\TP\external\glfw-3.1.2\tests && C:\msys64\ucrt64\bin\cc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S C:\Users\mateu\Documents\Bureau\prog_3d\Prog3D\Cameras_TP_21_11_2024\external\glfw-3.1.2\deps\getopt.c -o CMakeFiles\msaa.dir\__\deps\getopt.c.s

# Object files for target msaa
msaa_OBJECTS = \
"CMakeFiles/msaa.dir/msaa.c.obj" \
"CMakeFiles/msaa.dir/__/deps/getopt.c.obj"

# External object files for target msaa
msaa_EXTERNAL_OBJECTS =

external/glfw-3.1.2/tests/msaa.exe: external/glfw-3.1.2/tests/CMakeFiles/msaa.dir/msaa.c.obj
external/glfw-3.1.2/tests/msaa.exe: external/glfw-3.1.2/tests/CMakeFiles/msaa.dir/__/deps/getopt.c.obj
external/glfw-3.1.2/tests/msaa.exe: external/glfw-3.1.2/tests/CMakeFiles/msaa.dir/build.make
external/glfw-3.1.2/tests/msaa.exe: external/glfw-3.1.2/src/libglfw3.a
external/glfw-3.1.2/tests/msaa.exe: external/glfw-3.1.2/tests/CMakeFiles/msaa.dir/linkLibs.rsp
external/glfw-3.1.2/tests/msaa.exe: external/glfw-3.1.2/tests/CMakeFiles/msaa.dir/objects1.rsp
external/glfw-3.1.2/tests/msaa.exe: external/glfw-3.1.2/tests/CMakeFiles/msaa.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=C:\Users\mateu\Documents\Bureau\prog_3d\Prog3D\Cameras_TP_21_11_2024\TP\CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking C executable msaa.exe"
	cd /d C:\Users\mateu\Documents\Bureau\prog_3d\Prog3D\Cameras_TP_21_11_2024\TP\external\glfw-3.1.2\tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\msaa.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
external/glfw-3.1.2/tests/CMakeFiles/msaa.dir/build: external/glfw-3.1.2/tests/msaa.exe
.PHONY : external/glfw-3.1.2/tests/CMakeFiles/msaa.dir/build

external/glfw-3.1.2/tests/CMakeFiles/msaa.dir/clean:
	cd /d C:\Users\mateu\Documents\Bureau\prog_3d\Prog3D\Cameras_TP_21_11_2024\TP\external\glfw-3.1.2\tests && $(CMAKE_COMMAND) -P CMakeFiles\msaa.dir\cmake_clean.cmake
.PHONY : external/glfw-3.1.2/tests/CMakeFiles/msaa.dir/clean

external/glfw-3.1.2/tests/CMakeFiles/msaa.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" C:\Users\mateu\Documents\Bureau\prog_3d\Prog3D\Cameras_TP_21_11_2024 C:\Users\mateu\Documents\Bureau\prog_3d\Prog3D\Cameras_TP_21_11_2024\external\glfw-3.1.2\tests C:\Users\mateu\Documents\Bureau\prog_3d\Prog3D\Cameras_TP_21_11_2024\TP C:\Users\mateu\Documents\Bureau\prog_3d\Prog3D\Cameras_TP_21_11_2024\TP\external\glfw-3.1.2\tests C:\Users\mateu\Documents\Bureau\prog_3d\Prog3D\Cameras_TP_21_11_2024\TP\external\glfw-3.1.2\tests\CMakeFiles\msaa.dir\DependInfo.cmake "--color=$(COLOR)"
.PHONY : external/glfw-3.1.2/tests/CMakeFiles/msaa.dir/depend

