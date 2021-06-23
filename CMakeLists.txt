#
#	Created by WW on 2019-6-27, revised by WW on 2020-02-06.
#	Copyright © WW. All rights reserved.
#
#	This is a general template illustrates how to use Cmake.
#

cmake_minimum_required(VERSION 3.0)  # required cmake version at least

project(Phase_field) # project name

# Add directory "src" into the project, CMake will also execute the CMakeLists.txt file in "src",
# so you can use the files and libraries in "src".
add_subdirectory(src)
