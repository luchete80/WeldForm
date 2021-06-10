
project				(matrixmul LANGUAGES CXX CUDA)
ADD_EXECUTABLE  	(matrixmul "kernel.cu" "matrixmul.cu") 
	