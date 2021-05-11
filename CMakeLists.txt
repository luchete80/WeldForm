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
include_directories("D:/Luciano/Numerico/blitz/include")
link_directories("D:/Luciano/Numerico/blitz/lib")
  # using GCC
link_directories("D:/Luciano/Numerico/CompactNSearch_bin_mingw")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  # using Intel C++
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
MESSAGE("Using MSVC compilers")
# using Visual Studio C++
include_directories("D:/Luciano/Numerico/Libs/blitz-0.9_install_mswin/include")
include_directories("D:/Luciano/Numerico/Libs/gsl-2.5_install_mswin/include")
include_directories("D:/Luciano/Numerico/Libs/hdf5-1.12.0_install_mswin/include")
link_directories("D:/Luciano/Numerico/Libs/blitz-0.9_install_mswin/lib")
link_directories("D:/Luciano/Numerico/CompactNSearch_bin_nmake_release")

include_directories("Source")
include_directories("External")

add_subdirectory("examples")
endif()

include_directories("D:/Luciano/Numerico/CompactNSearch/include")


SET(EXES
	# 1-Poiseuille
	# 2-Couette
	# 3-DamBreak
	# 4-DamBreak3D
	#Cylinder_Test
	#Plate-Yield
	Compression
	#Plate-Yield_hcorr
	#Plate-Yield_CNS
	#Plate-Yield_Hard
	#Plate-Yield_CNS_Hard
	#Plate-Yield_CNS_v2
	#Plate-Yield_3d_al
	#Plate-Yield_3d_al_CNS
	#SolidBlock_Heat
	#SolidBlock_Heat_CNS
	#Neighbour_Test
	#Neighbour_Test_CNS
)


FOREACH(var ${EXES})
    ADD_EXECUTABLE        (${var} "${var}.cpp" 
	"${CMAKE_SOURCE_DIR}/Source/CompactNSearch/CompactNSearch.cpp"
	)
    TARGET_LINK_LIBRARIES (${var} ${LIBS} blitz)
    SET_TARGET_PROPERTIES (${var} PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")
ENDFOREACH(var)

    # ADD_EXECUTABLE        (${EXES} "Plate-Yield.cpp")
    # TARGET_LINK_LIBRARIES (${EXES} ${LIBS} blitz)
    # SET_TARGET_PROPERTIES (${EXES} PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")
