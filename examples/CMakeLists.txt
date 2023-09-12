add_subdirectory(input)


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
  #Contact_Compression_2d
  #test_new_sym
  #Cantilever
  # examples/Compression
	#Compression_XSPH
  #Compression_AVisc
  #Compression_1010
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
  #Compression_thermal
  #Compression_nastran
  #Contact_Compression_thermal
  #Contact_Compression_Rot
	#Cantilever
	#Cantilever2
	#Cantilever2_pl
	#Contact_Neighbours
	# examples/Contact_Compression
  # examples/Contact_Compression_1010
  #Contact_Compression_1010_1side
  #Contact_Compression_fast_short
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
  #Contact_Compression_1010_unif
  #Compression_1010
  # Contact_Compression_1010
  #Compression_1010_qt
  #Contact_Compression_1010_half
  #Compression_1010_half_ghost
	#NastranReader_test
	#CompressionForces
  # Compression 
  #Compression_2d
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
  #FSW_red
  #Drilling
  #FSW_nastran
	#Elastic
	#Heat_2D
	#Impact_Campos
	#Impact_Contact
	Impact_Zhang2017
  #Kernel_Test
	#Kernel_Test_2D
	#Kernel_Test_2D_vel
	#Traction
	#traction_app
	#Traction_Hollomon
	#Traction_JC
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
  #SolidBlock_Heat_T
	#slice_test
  #SolidBlock_HeatSOA
	#SolidBlock_Heat_CNS
	#SolidBlock_Heat_conv3faces
	#SolidBlock_Heat_source
	#SolidBlock_Heat_conv3faces_Flux
	#SolidBlock_Heat_source2
  #SolidBlock_Heat_T
	#SolidBlock_Heat_source2_boundpart
	#Thermal_Expansion_2
  #Thick_Cyl_AxiSymm
	
	#Neighbour_Test
	#Neighbour_Test_CNS
	
  #"restart_test"
)

FOREACH(var ${EXES})
    ADD_EXECUTABLE        (${var} "${var}.cpp" 
	#"${CMAKE_SOURCE_DIR}/Source/Mesh.cpp"	#Contact Problems
  #"${CMAKE_SOURCE_DIR}/Source/NastranVolReader.cpp"
  
	)
    TARGET_LINK_LIBRARIES (${var} ${LIBS} blitz)
    SET_TARGET_PROPERTIES (${var} PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")
ENDFOREACH(var)