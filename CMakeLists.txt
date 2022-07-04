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
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  # using Intel C++
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
MESSAGE("Using MSVC compilers")
# using Visual Studio C++
include_directories("D:/Luciano/Numerico/Libs/blitz-0.9_install_mswin/include")
include_directories("D:/Luciano/Numerico/Libs/gsl-2.5_install_mswin/include")
include_directories("D:/Luciano/Numerico/Libs/hdf5-1.12.0_install_mswin/include")
link_directories("D:/Luciano/Numerico/Libs/blitz-0.9_install_mswin/lib")

include_directories("Source")
include_directories("External")

endif()

if (GPU_NBSEARCH)
add_definitions(-D_GPU_NEIBSEARCH)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GPU_NEIBSEARCH")

endif()

#add_subdirectory("examples")

include_directories("D:/Luciano/Numerico/CompactNSearch/include")


SET(EXES
	#Ale
	#Ale_modelo6
	#SolidBlock_Heat_conv3faces_Flux_n30
	#Ale_Modelo6_placa
	# 1-Poiseuille
	# 2-Couette
	# 3-DamBreak
	# 4-DamBreak3D
	#Cylinder_Test
	#Plate-Yield
	#Compression_XSPH
  #Compression_AVisc
  #Compression_half_ghost_vel
	#Compression_thermal
	#Compression_12_solver_orig
	#Compression_smallD
	#Compression_smallL
	#Compression_smallL_1side
	#Compression_smallD_CNS
	#Compression_10ts
	#Compression_Elastic_CNS_10ts
	#Compression_Elastic_CNS
  #Compression_TRIP_half_ghost_vel
  #Compression_TRIP_vel
	#Cantilever
	#Cantilever2
	#Cantilever2_pl
	#Contact_Neighbours
	#Contact_Compression
  #Contact_Compression_2surf
  #Contact_Compression_vel
  #Contact_Compression_test
  #Contact_Compression_half_ghost
  #Contact_Compression_half_ghost_debug
  #Contact_Compression_qt_ghost
  #Contact_Compression_TRIP
  #Contact_Compression_TRIP_half_ghost
  #Contact_Compression_TRIP_2surf
  #Contact_Compression_1010_hlaf_ghost
  #Compression_1010
  #Contact_Compression_1010
  #Compression_1010_half_ghost
	#NastranReader_test
	#CompressionForces
  #Compression
	#CompressionVel
  #CompressionVel_Implicit
	#CompressionForces_Ghost
	#Compression_Free
  #Compression_Free_GPU_v1
  #Compression_qt
	#Compression_qt_ghost
  #Compression_half  
  #Compression_half_ghost_vel
  #CompressionVel_2sides
  #Compression_half_ghost
	#CompressionForces2
	#Cantilever_comp_contact
	#Contact_Cantilever
	#FreeSurface
  #FSW
	#Elastic
	#Heat_2D
	#Impact_Campos
	#Impact_Contact
	#Impact_Zhang2017
  #Kernel_Test
	#Kernel_Test_2D
	#Kernel_Test_2D_vel
	#Traction
	#traction_app
	#Traction_Hollomon
	#Traction_JC
	#ThermalCompression
	#Thermal_Expansion_2
	#4-DamBreak3D
	#Plate-Yield
	#Plate-Yield_hcorr
	#Plate-Yield_CNS
	#Plate-Yield_Hard
	#Plate-Yield_CNS_Hard
	#Plate-Yield_CNS_v2
	#Plate-Yield_3d_al
	#Plate-Yield_3d_al_CNS
	#SolidBlock_HeatSOA
	#SolidBlock_Heat_CNS
	#SolidBlock_Heat_conv3faces
	#SolidBlock_Heat_source
	#SolidBlock_Heat_conv3faces_Flux
	#SolidBlock_Heat_source2
	#SolidBlock_Heat_source2_boundpart
	#Thermal_Expansion_2
	
	#Neighbour_Test
	#Neighbour_Test_CNS
	#"Source/WeldForm"
)

add_subdirectory("examples")

FOREACH(var ${EXES})
    ADD_EXECUTABLE        (${var} "${var}.cpp" 
	#"${CMAKE_SOURCE_DIR}/Source/Mesh.cpp"	#Contact Problems
	)
    TARGET_LINK_LIBRARIES (${var} ${LIBS} blitz)
    SET_TARGET_PROPERTIES (${var} PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")
ENDFOREACH(var)

ADD_EXECUTABLE        (WeldForm "${PROJECT_SOURCE_DIR}/Source/WeldForm.cpp" )
TARGET_LINK_LIBRARIES (WeldForm ${LIBS} blitz)
SET_TARGET_PROPERTIES (WeldForm PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")	
    # ADD_EXECUTABLE        (${EXES} "Plate-Yield.cpp")
    # TARGET_LINK_LIBRARIES (${EXES} ${LIBS} blitz)
    # SET_TARGET_PROPERTIES (${EXES} PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")