#
#	Created by WW on 2019-6-22.
#	Copyright © WW. All rights reserved.
#
#	This is a general template illustrates how to use Cmake.
#

cmake_minimum_required(VERSION 3.0)  # required cmake version at least

# project(Phase_field) # project name

# Bring the header-files in directory "include" into the project, 
# which equals to "-I<inc>" flag of gcc,
# so you don't need to include these header-files in the "add_executable/library" Command.
include_directories(/home/wwang138/.local/include)

# Link current CMakeLists.txt directory, 
# then you can target_link_libraries() in these directories.
# You can also link other directories, which equals to "-L<path>" flag of gcc.

# Add directory "src" into the project, CMake will also execute the CMakeLists.txt file in "src",
# so you can use the files and libraries in "src".
# add_subdirectory(src)

# Cmake C/C++ compiler will use C++ 11, which equals to "set(CMAKE_CXX_STANDARD 11)".	
add_compile_options(-std=c++11)

# Add some new flags into original CXXFLAGS.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")


# Define some variables representing the files.
# GLOB and GLOB_RECURSE can recognize the Regular Expression "*".
# file(GLOB LIB_HEADERS *.h *.hpp)
# file(GLOB LIB_SOURCES *.cpp *.c)

# Compile and creat new library (SHARED or STATIC).
# add_library(MYLIB STATIC ${LIB_HEADERS} ${LIB_SOURCES})

# set(EXTRA_LIBS ${EXTRA_LIBS} MYLIB)  # add libraries into the EXTRA_LIBS
find_package(MPI REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C CXX)

include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${MPI_INCLUDE_DIRS})
# Define some variables representing the files.
# GLOB and GLOB_RECURSE can recognize the Regular Expression "*".
file(GLOB HEADERS *.h *.hpp)
file(GLOB SOURCES *.c *.cpp)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR})
add_executable(test.exe ${SOURCES} ${HEADERS})

# Link the executable file with all libraries
# target_link_libraries(magnon.exe ${EXTRA_LIBS})

# If you already have some libraries, you could link them after link_directories(),
# which equals to -l<libname>
target_link_libraries(test.exe /home/wwang138/.local/lib/libnrutil.a)
target_link_libraries(test.exe ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})
target_link_libraries(test.exe ${MPI_CXX_LIBRARIES} ${MPI_LIBRARIES})
