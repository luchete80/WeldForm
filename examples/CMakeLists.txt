project(cuda_sum LANGUAGES CXX CUDA)
# Visual studio solution directories.
set_property(GLOBAL PROPERTY USE_FOLDERS on)

add_executable(cuda_sum cuda_sum.cu)
add_executable(saxpy saxpy.cu)
target_include_directories(cuda_sum PUBLIC
	${CUDA_INCLUDE_DIRS}
)
target_include_directories(saxpy PUBLIC
	${CUDA_INCLUDE_DIRS}
)
target_link_libraries(cuda_sum PUBLIC ${CUDA_LIBRARIES})
target_link_libraries(saxpy PUBLIC ${CUDA_LIBRARIES})