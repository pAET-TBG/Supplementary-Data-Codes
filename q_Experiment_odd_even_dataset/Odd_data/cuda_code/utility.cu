#include <iostream>
#include "mex.h"
#include "omp.h"
#include <cmath>
#include "gpu/mxGPUArray.h"
#include <cuda_runtime.h>
using namespace std;

#define FLOAT_MAX 1e5

static const int blockSize = 1024;
static const int gridSize = 24; 

template <typename T>
__global__ void maxCommMultiBlock(const T *gArr, long long arraySize, T *gOut) {
    int thIdx = threadIdx.x;
    int gthIdx = thIdx + blockIdx.x*blockSize;
    const int gridSize = blockSize*gridDim.x;
    T max_val = -FLOAT_MAX;
    for (int i = gthIdx; i < arraySize; i += gridSize)
        max_val = max(max_val, gArr[i] );
    __shared__ T shArr[blockSize];
    shArr[thIdx] = max_val;
    __syncthreads();
    for (int size = blockSize/2; size>0; size/=2) { //uniform
        if (thIdx<size)
            shArr[thIdx] = max( shArr[thIdx], shArr[thIdx+size] );
        __syncthreads();
    }
    if (thIdx == 0){
        gOut[blockIdx.x] = shArr[0];
        //atomicMaxFloat(gOut, shArr[0]);
    }
}

template <typename T>
__host__ double maxArray(T* arr_pt, long long wholeArraySize) {
    //T* arr_pt;
    //cudaMalloc((void**)&arr_pt, wholeArraySize * sizeof(T));
    //cudaMemcpy(arr_pt, arr, wholeArraySize * sizeof(T), cudaMemcpyHostToDevice);

    T sum;
    T* sum_pt;
    cudaMalloc((void**)&sum_pt, sizeof(T)*gridSize);
    
    maxCommMultiBlock<<<gridSize, blockSize>>>(arr_pt, wholeArraySize, sum_pt);
    //dev_out now holds the partial result
    maxCommMultiBlock<<<1, blockSize>>>(sum_pt, gridSize, sum_pt);
    //dev_out[0] now holds the final result
    cudaDeviceSynchronize();
    
    cudaMemcpy(&sum, sum_pt, sizeof(T), cudaMemcpyDeviceToHost);
    //cudaFree(arr_pt);
    cudaFree(sum_pt);
    return sum;
}



template <typename T>
__global__ void minCommMultiBlock(const T *gArr, long long arraySize, T *gOut) {
    int thIdx = threadIdx.x;
    int gthIdx = thIdx + blockIdx.x*blockSize;
    const int gridSize = blockSize*gridDim.x;
    T max_val = FLOAT_MAX;
    for (int i = gthIdx; i < arraySize; i += gridSize)
        max_val = min(max_val, gArr[i] );
    __shared__ T shArr[blockSize];
    shArr[thIdx] = max_val;
    __syncthreads();
    for (int size = blockSize/2; size>0; size/=2) { //uniform
        if (thIdx<size)
            shArr[thIdx] = min( shArr[thIdx], shArr[thIdx+size] );
        __syncthreads();
    }
    if (thIdx == 0){
        gOut[blockIdx.x] = shArr[0];
        //atomicMinFloat(gOut, shArr[0]);
    }
}

template <typename T>
__host__ double minArray(T* arr_pt, long long wholeArraySize) {
    //T* arr_pt;
    //cudaMalloc((void**)&arr_pt, wholeArraySize * sizeof(T));
    //cudaMemcpy(arr_pt, arr, wholeArraySize * sizeof(T), cudaMemcpyHostToDevice);

    T sum;
    T* sum_pt;
    cudaMalloc((void**)&sum_pt, sizeof(T)*gridSize);
    
    minCommMultiBlock<<<gridSize, blockSize>>>(arr_pt, wholeArraySize, sum_pt);
    //dev_out now holds the partial result
    minCommMultiBlock<<<1, blockSize>>>(sum_pt, gridSize, sum_pt);
    //dev_out[0] now holds the final result
    cudaDeviceSynchronize();
    
    cudaMemcpy(&sum, sum_pt, sizeof(T), cudaMemcpyDeviceToHost);
    //cudaFree(arr_pt);
    cudaFree(sum_pt);
    return sum;
}


template <typename T>
__global__ void sumCommMultiBlock(const T *gArr, long long arraySize, T *gOut) {
    int thIdx = threadIdx.x;
    int gthIdx = thIdx + blockIdx.x*blockSize;
    const int gridSize = blockSize*gridDim.x;
    T sum = 0;
    for (int i = gthIdx; i < arraySize; i += gridSize)
        sum += abs(gArr[i]);
    __shared__ T shArr[blockSize];
    shArr[thIdx] = sum;
    __syncthreads();
    for (int size = blockSize/2; size>0; size/=2) { //uniform
        if (thIdx<size)
            shArr[thIdx] += shArr[thIdx+size];
        __syncthreads();
    }
    if (thIdx == 0)
        gOut[blockIdx.x] = shArr[0];
}

template <typename T>
__host__ double sumArray(const T* arr, long long wholeArraySize) {
    T* arr_pt;
    cudaMalloc((void**)&arr_pt, wholeArraySize * sizeof(T));
    cudaMemcpy(arr_pt, arr, wholeArraySize * sizeof(T), cudaMemcpyHostToDevice);

    T sum;
    T* sum_pt;
    cudaMalloc((void**)&sum_pt, sizeof(T)*gridSize);
    
    sumCommMultiBlock<<<gridSize, blockSize>>>(arr_pt, wholeArraySize, sum_pt);
    //dev_out now holds the partial result
    sumCommMultiBlock<<<1, blockSize>>>(sum_pt, gridSize, sum_pt);
    //dev_out[0] now holds the final result
    cudaDeviceSynchronize();
    
    cudaMemcpy(&sum, sum_pt, sizeof(T), cudaMemcpyDeviceToHost);
    cudaFree(arr_pt);
    cudaFree(sum_pt);
    return sum;
}