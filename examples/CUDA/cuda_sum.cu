#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <iostream>
#include <numeric>
using namespace std;

__global__ void sum(int* input)
{
	const int tid = threadIdx.x;

	auto step_size = 1;
	int number_of_threads = blockDim.x;
	
	printf("Hello from block %d, thread %d\n", blockIdx.x, threadIdx.x);
	while (number_of_threads > 0)
	{
		if (tid < number_of_threads) // still alive?
		{
			const auto fst = tid * step_size * 2;
			const auto snd = fst + step_size;
			input[fst] += input[snd];
		}

		step_size <<= 1; 
		number_of_threads >>= 1;
	}
}

int main()
{
	const auto count = 16;
	const int size = count * sizeof(int);
	//int h[] = {13, 27, 15, 14, 33, 2, 24, 6};
	
	
	int h[] = {13, 27, 15, 14, 33, 2, 24, 6, 1, 120, 40, 50 , 60, 12, 13, 18};
	
	int* d;
	
	cudaMalloc(&d, size);
	cudaMemcpy(d, h, size, cudaMemcpyHostToDevice);

	sum <<<1, count / 2 >>>(d);

	int result;
	cudaMemcpy(&result, d, sizeof(int), cudaMemcpyDeviceToHost);

	cout << "Sum is " << result << endl;

	getchar();

	cudaFree(d);
	delete[] h;

	return 0;
}