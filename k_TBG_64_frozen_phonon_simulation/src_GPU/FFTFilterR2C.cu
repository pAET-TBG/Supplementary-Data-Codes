#include <iostream>
#include <stdio.h>
#include "mex.h"
#include "omp.h"
#include <cmath>
#include "gpu/mxGPUArray.h"
#include <cuda_runtime.h>
#include <cufft.h>

// #ifndef CUDA_CHECK_DEFINED
// #define CUDA_CHECK_DEFINED
// #define CUDA_CHECK(err) if(err!=cudaSuccess){mexErrMsgIdAndTxt("CUDA",cudaGetErrorString(err));}
// #endif

#ifndef CUDA_CHECK_DEFINED
#define CUDA_CHECK_DEFINED
#define CUDA_CHECK(call) do { \
    cudaError_t err = (call); \
    if (err != cudaSuccess) { \
        mexErrMsgIdAndTxt("cudaError", cudaGetErrorString(err)); \
    } \
} while(0)
#endif // CUDA_CHECK_DEFINED

#define CUFFT_CHECK(err) if(err!=CUFFT_SUCCESS){mexErrMsgIdAndTxt("CUFFT","cuFFT error code %d",(int)err);}

template<typename T>
__global__ void normalizeKernel(const T* in, T* out, const size_t N, const size_t scale) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid < N) {
        out[tid] = in[tid] / scale;
    }
}

// Promote real to complex (imag = 0)
__global__ void realToComplexKernel(const float* __restrict__ in,
                                    cufftComplex* __restrict__ out,
                                    size_t total)
{
    size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid < total) {
        out[tid].x = in[tid];
        out[tid].y = 0.0f;
    }
}


__global__ void get_B_half(const float* __restrict__ B,
                                float* __restrict__ B_half,
                        int N1, int N2, int N1C)
{
    int ind1 = blockIdx.x * blockDim.x + threadIdx.x;
    int ind2 = blockIdx.y * blockDim.y + threadIdx.y;

    if (ind1<N1C && ind2<N2){
        B_half[ind1 + N1C*ind2] = B[ind1 + N1*ind2];
    }
}    

template <typename T>
__global__ void circShift2D(const T* __restrict__ data_in,
                            T* __restrict__ data_out,
                            int N1, int N2,
                            int num_slices,
                            int sh_1, int sh_2)
{
    // MATLAB-style row/column indices
    int ind_1 = blockIdx.x * blockDim.x + threadIdx.x;   // row index
    int ind_2 = blockIdx.y * blockDim.y + threadIdx.y;   // col index

    if (ind_1 >= N1 || ind_2 >= N2)
        return;

    int new_1 = (ind_1+sh_1+N1) % N1;
    int new_2 = (ind_2+sh_2+N2) % N2;

    // Base pointers for slice iteration
    int idx_in  = ind_1 + ind_2 * N1;
    int idx_out   = new_1 + new_2 * N1;
    int sliceSize = N1 * N2;

    // Copy this (x,y) pixel across all slices
    for (int k = 0; k < num_slices; k++) {
        data_out[k*sliceSize + idx_out] =
        data_in [k*sliceSize + idx_in];
    }
}    

// Complex * real (B is 2D, broadcast over num_pj)
__global__ void multiplyComplexReal(const cufftComplex* __restrict__ A,
                                    const float* __restrict__ B,
                                    cufftComplex* __restrict__ C,
                                    int N1, int N2, int num_pj)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int sliceSize = N1 * N2;
    int total = sliceSize * num_pj;

    if (tid < total) {
        int idx2D = tid % sliceSize;   // broadcast B over projections
        float b   = B[idx2D];
        cufftComplex a = A[tid];
        C[tid].x = a.x * b;
        C[tid].y = a.y * b;
    }
}

class FFTFilterR2C
{
public:
    using Real    = float;
    using Complex = cufftComplex;

    FFTFilterR2C(const mxArray* B_in, int num_pj)
        : FFTFilterR2C( mxGPUCreateFromMxArray(B_in) , num_pj, true)
    {}

    FFTFilterR2C(const mxGPUArray* B_gpu, int num_pj)
        : FFTFilterR2C( B_gpu , num_pj, false)
    {}

    FFTFilterR2C(const mxGPUArray* B_gpu, int num_pj, bool createNewGPU)
        : num_pj_(num_pj),
          N1_(0), N2_(0), N1c_(0),
          sliceSize_(0), sliceSizeCpl_(0),
          totalRealA_(0), totalCplxA_(0),
          dB_shift_(nullptr), dB_half_(nullptr),
          dA_tmp_(nullptr), dA_fft_(nullptr),
          plan_r2c_(0), plan_c2r_(0)
    {
        // -----------------------------
        // Validate B and derive sizes
        // -----------------------------
        if (mxGPUGetNumberOfDimensions(B_gpu) != 2)
            mexErrMsgIdAndTxt("invalidB","B must be [N1,N2]");

        const mwSize* dimsB = mxGPUGetDimensions(B_gpu);
        N1_ = static_cast<int>(dimsB[0]);
        N2_ = static_cast<int>(dimsB[1]);
        N1c_ = N1_/2 + 1;

        sliceSize_    = static_cast<size_t>(N1_) * N2_;
        sliceSizeCpl_ = static_cast<size_t>(N2_) * N1c_;
        totalRealA_   = sliceSize_ * num_pj_;
        totalCplxA_   = sliceSizeCpl_ * num_pj_;

        const Real* dB_in = static_cast<const Real*>(mxGPUGetDataReadOnly(B_gpu));
        if(createNewGPU){
            mxGPUDestroyGPUArray(B_gpu);
        }        

        // const Real* hB_in = static_cast<const Real*>(mxGetData(B_in));
        // Real* dB_in;
        // cudaMalloc(&dB_in, sliceSize_*sizeof(Real) );
        // cudaMemcpy(dB_in, hB_in, sliceSize_*sizeof(Real), cudaMemcpyHostToDevice);
        // CUDA_CHECK(cudaGetLastError());

        CUDA_CHECK(cudaMalloc(&dB_shift_, sliceSize_    * sizeof(Real)));
        CUDA_CHECK(cudaMalloc(&dB_half_,  sliceSizeCpl_ * sizeof(Real)));        

        // -----------------------------
        // Precompute ifftshift(B)
        // -----------------------------
        sh_ifft1_ = (N1_ + 1) / 2;
        sh_ifft2_ = (N2_ + 1) / 2;
        sh_fft1_  =  N1_ / 2;
        sh_fft2_  =  N2_ / 2;

        // dim3 block2D(16,16);
        // dim3 grid2D((N1_+15)/16, (N2_+15)/16);
        block2D = dim3(16,16);
        grid2D  = dim3((N1_+15)/16, (N2_+15)/16);        
        threads = 256;
        blocks_cplxA  = (totalCplxA_ + threads - 1) / threads;  
        blocks_realA  = (totalRealA_ + threads - 1) / threads;   
        
        circShift2D<Real><<<grid2D,block2D>>>(
            dB_in, dB_shift_,
            N1_, N2_, 1,
            sh_ifft1_, sh_ifft2_);
        CUDA_CHECK(cudaGetLastError());
        dB_in = nullptr;
        
        // CUDA_CHECK(cudaDeviceSynchronize());  // REQUIRED
        // cudaFree(dB_in);
        // CUDA_CHECK(cudaGetLastError());

        // -----------------------------
        // Extract half-spectrum of B
        // -----------------------------
        dim3 blockB(16,16);
        dim3 gridB((N1c_+15)/16, (N2_+15)/16);

        get_B_half<<<gridB,blockB>>>(
            dB_shift_, dB_half_, N1_, N2_, N1c_);
        CUDA_CHECK(cudaGetLastError());

        // -----------------------------
        // Allocate reusable A buffers
        // -----------------------------
        CUDA_CHECK(cudaMalloc(&dA_tmp_, totalRealA_ * sizeof(Real)));
        CUDA_CHECK(cudaMalloc(&dA_fft_, totalCplxA_ * sizeof(Complex)));

        // -----------------------------
        // Create FFT plans
        // -----------------------------
        int n[2] = { N2_, N1_ };  // slow, fast

        int inembed[2] = { N2_, N1_ };
        int onembed[2] = { N2_, N1c_ };

        size_t istride = 1;
        size_t ostride = 1;
        size_t idist   = sliceSize_;
        size_t odist   = sliceSizeCpl_;

        CUFFT_CHECK(cufftPlanMany(&plan_r2c_, 2, n,
                                  inembed, istride, idist,
                                  onembed, ostride, odist,
                                  CUFFT_R2C, num_pj_));

        CUFFT_CHECK(cufftPlanMany(&plan_c2r_, 2, n,
                                  onembed, ostride, odist,
                                  inembed, istride, idist,
                                  CUFFT_C2R, num_pj_));
        // mxGPUDestroyGPUArray(B_gpu);
    }

    ~FFTFilterR2C()
    {
        if (plan_r2c_) CUFFT_CHECK(cufftDestroy(plan_r2c_));
        if (plan_c2r_) CUFFT_CHECK(cufftDestroy(plan_c2r_));
        
        if (dB_shift_) CUDA_CHECK(cudaFree(dB_shift_));
        if (dB_half_)  CUDA_CHECK(cudaFree(dB_half_));
        if (dA_tmp_)   CUDA_CHECK(cudaFree(dA_tmp_));
        if (dA_fft_)   CUDA_CHECK(cudaFree(dA_fft_));
    }

    // ------------------------------------------------------------
    // apply(): compute
    // C = fftshift( ifft2( fft2(ifftshift(A)) .* ifftshift(B) ) )
    // ------------------------------------------------------------
    void apply(const Real* dA_in, Real* dC_out)
    {
        // dim3 block2D(16,16);
        // dim3 grid2D((N1_+15)/16, (N2_+15)/16);
        // -----------------------------
        // 1. ifftshift(A)
        // -----------------------------
        circShift2D<Real><<<grid2D,block2D>>>(
            dA_in, dA_tmp_,
            N1_, N2_, num_pj_,
            sh_ifft1_, sh_ifft2_);
        CUDA_CHECK(cudaGetLastError());


        // -----------------------------
        // 2. Forward FFT
        // -----------------------------
        CUFFT_CHECK(cufftExecR2C(
            plan_r2c_,
            reinterpret_cast<cufftReal*>(dA_tmp_),
            reinterpret_cast<cufftComplex*>(dA_fft_)));


        // -----------------------------
        // 3. Multiply F .* B_half
        // -----------------------------
        // int threads = 256;
        // int blocks  = (totalCplxA_ + threads - 1) / threads;
        multiplyComplexReal<<<blocks_cplxA,threads>>>(
            dA_fft_,
            dB_half_,
            dA_fft_,
            N1c_, N2_, num_pj_);
        CUDA_CHECK(cudaGetLastError());


        // -----------------------------
        // 4. Inverse FFT
        // -----------------------------
        CUFFT_CHECK(cufftExecC2R(
            plan_c2r_,
            reinterpret_cast<cufftComplex*>(dA_fft_),
            reinterpret_cast<cufftReal*>(dA_tmp_)));


        // -----------------------------
        // 5. fftshift
        // -----------------------------
        circShift2D<Real><<<grid2D,block2D>>>(
            dA_tmp_, dC_out,
            N1_, N2_, num_pj_,
            sh_fft1_, sh_fft2_);
        CUDA_CHECK(cudaGetLastError());

        
        // -----------------------------
        // 6. Normalize
        // -----------------------------
        // int blocksA = (totalRealA_ + threads - 1) / threads;
        normalizeKernel<<<blocks_realA,threads>>>(
            dC_out, dC_out, totalRealA_, sliceSize_);
        CUDA_CHECK(cudaGetLastError());
    }

    int N1() const { return N1_; }
    int N2() const { return N2_; }
    int numProj() const { return num_pj_; }
        

private:
    int num_pj_;
    int N1_, N2_, N1c_;
    size_t sliceSize_, sliceSizeCpl_;
    size_t totalRealA_, totalCplxA_;

    int sh_ifft1_, sh_ifft2_;
    int sh_fft1_,  sh_fft2_;

    Real*    dB_shift_;
    Real*    dB_half_;

    Real*    dA_tmp_;
    Complex* dA_fft_;

    dim3 block2D;
    dim3 grid2D;
    
    int threads;
    int blocks_cplxA;
    int blocks_realA;

    cufftHandle plan_r2c_;
    cufftHandle plan_c2r_;
    
};