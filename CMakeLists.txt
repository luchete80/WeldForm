#####################################################################################
# PersianSPH - A C++ library to simulate Mechanical Systems (solids, fluids         #
#             and soils) using Smoothed Particle Hydrodynamics method               #
# Copyright (C) 2016 Maziar Gholami Korzani and Sergio Galindo-Torres               #
#                                                                                   #
# This file is part of PersianSPH                                                   #
#                                                                                   #
# This is free software; you can redistribute it and/or modify it under the         #
# terms of the GNU General Public License as published by the Free Software         #
# Foundation; either version 3 of the License, or (at your option) any later        #
# version.                                                                          #
#                                                                                   #
# This program is distributed in the hope that it will be useful, but WITHOUT ANY   #
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A   #
# PARTICULAR PURPOSE. See the GNU General Public License for more details.          #
#                                                                                   #
# You should have received a copy of the GNU General Public License along with      #
# PersianSPH; if not, see <http://www.gnu.org/licenses/>                            #
#####################################################################################

CMAKE_MINIMUM_REQUIRED (VERSION 2.8)

PROJECT                (New_Test)
cmake_policy(SET CMP0054 NEW) #elseif
add_definitions(-D_CRT_NO_VA_START_VALIDATION)

add_definitions(-D_USE_MATH_DEFINES) #COMPACTNSEARCH??

INCLUDE ($ENV{SPH}/Modules/FindPKG.cmake)

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  # using Clang
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
MESSAGE("Using GNU compiler")
include_directories("$ENV{PKG}/blitz/include")
include_directories("$ENV{PKG}/blitz-0.9/include")
link_directories("$ENV{PKG}/blitz/lib")
link_directories("$ENV{PKG}/blitz-0.9/lib")
  # using GCC
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  # using Intel C++
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
MESSAGE("Using MSVC compilers")
# using Visual Studio C++
include_directories("$ENV{PKG}/blitz-0.9_install_mswin/include")
include_directories("$ENV{PKG}/gsl-2.5_install_mswin/include")
include_directories("$ENV{PKG}/hdf5-1.12.0_install_mswin/include")
link_directories("$ENV{PKG}/blitz-0.9_install_mswin/lib")

include_directories("Source")
include_directories("External")

endif()

if (GPU_NBSEARCH)
add_definitions(-D_GPU_NEIBSEARCH)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GPU_NEIBSEARCH")

endif()

add_subdirectory("examples")

include_directories("D:/Luciano/Numerico/CompactNSearch/include")


ADD_EXECUTABLE        (WeldForm "Source/WeldForm.cpp" 
)

add_subdirectory ("${PROJECT_SOURCE_DIR}/External/LSDynaReader/src")

TARGET_LINK_LIBRARIES (WeldForm ${LIBS} blitz lsdynareader)
SET_TARGET_PROPERTIES (WeldForm PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")



# ADD_EXECUTABLE        (WeldForm "${PROJECT_SOURCE_DIR}/Source/WeldForm.cpp" )
# TARGET_LINK_LIBRARIES (WeldForm ${LIBS} blitz)
# SET_TARGET_PROPERTIES (WeldForm PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")	


