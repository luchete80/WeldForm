project(cuda_sum LANGUAGES CXX CUDA)
# Visual studio solution directories.
set_property(GLOBAL PROPERTY USE_FOLDERS on)

add_executable(cuda_sum cuda_sum.cu)
target_include_directories(cuda_sum PUBLIC
	${CUDA_INCLUDE_DIRS}
)
target_link_libraries(cuda_sum PUBLIC ${CUDA_LIBRARIES})
