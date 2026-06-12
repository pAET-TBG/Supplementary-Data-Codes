#include <iostream>
#include <stdio.h>
#include "mex.h"
#include "omp.h"
#include <cmath>
#include "gpu/mxGPUArray.h"
// #include "FFTFilterR2C.cu"
#include "TomoModelOptimizer.cu"
using namespace std;

#ifdef printf
#undef printf
#endif

// #ifndef CUDA_CHECK
// #define CUDA_CHECK() { \
//     cudaError_t err = cudaGetLastError(); \
//     if (err != cudaSuccess) { \
//         mexErrMsgIdAndTxt("accumulateProjection:cudaError", cudaGetErrorString(err)); \
//     } \
// }
// #endif

#ifndef CUDA_CHECK_DEFINED
#define CUDA_CHECK_DEFINED
#define CUDA_CHECK(call) do { \
    cudaError_t err = (call); \
    if (err != cudaSuccess) { \
        mexErrMsgIdAndTxt("accumulateProjection:cudaError", cudaGetErrorString(err)); \
    } \
} while(0)
#endif


#if __CUDA_ARCH__ < 600
template <typename T>
__device__ double atomicAdd(T* address, T val)
{
    unsigned long long int* address_as_ull =
    (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
        __double_as_longlong(val +
        __longlong_as_double(assumed)));

        // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);

    return __longlong_as_double(old);
}
#endif


__device__ float roundAwayFromZero(float x) {
    if (x >= 0.0f) {
        return floorf(x + 0.5f);
    } else {
        return ceilf(x - 0.5f);
    }
}

__device__ float matlabRound(float x) {
    float r = truncf(x);              // truncate toward zero
    float diff = fabsf(x - r);        // fractional part
    if (diff >= 0.5f) {
        // move one step away from zero
        r += copysignf(1.0f, x);
    }
    return r;
}


// __global__ void forward_kernel(
//     const float* model,
//     const float* Rs,
//     const int* atoms, int num_atoms,
//     const float* b, const float* h,
//     float* proj,
//     // float* Grad_h,
//     // float* Grad_b,
//     float* Grad_x,
//     float* Grad_y,
//     float* Grad_z,
//     int N1, int N2, int halfWidth,
//     int num_pj)
// {
//     int jj = blockIdx.x;  
//     int ii = blockIdx.y * blockDim.y + threadIdx.y;  

//     if (jj >= num_pj || ii >= num_atoms) return;


//     const float* Rs_jj = Rs + 9*jj;

//     float X_cen = Rs_jj[0]*model[3*ii] + Rs_jj[3]*model[3*ii+1] + Rs_jj[6]*model[3*ii+2];
//     float Y_cen = Rs_jj[1]*model[3*ii] + Rs_jj[4]*model[3*ii+1] + Rs_jj[7]*model[3*ii+2];
//     float Z_cen = Rs_jj[2]*model[3*ii] + Rs_jj[5]*model[3*ii+1] + Rs_jj[8]*model[3*ii+2];

//     // int X_round = __float2int_rn(X_cen);
//     // int Y_round = __float2int_rn(Y_cen);
//     // int Z_round = __float2int_rn(Z_cen);

//     float X_round = roundf(X_cen);
//     float Y_round = roundf(Y_cen);
//     float Z_round = roundf(Z_cen);    

//     // float X_round = roundAwayFromZero(X_cen);
//     // float Y_round = roundAwayFromZero(Y_cen);
//     // float Z_round = roundAwayFromZero(Z_cen);     
    
//     int X_round_int = X_round;
//     int Y_round_int = Y_round;
//     // int Z_round_int = static_cast<int>(Z_round);

//     // if(jj==1 && ii == 9841){
//     //     printf("%d,%d: %f,%f,%f\n",jj,ii, X_cen, Y_cen, Z_cen);
//     //     printf("%d,%d: %f,%f,%f\n",jj,ii, X_round, Y_round, Z_round);
//     //     printf("%d,%d: %d,%d\n",jj,ii, X_round_int, Y_round_int);
//     // }

//     int N1_half = N1/2;
//     int N2_half = N2/2;

//     int b_size = 2*halfWidth+1;

//     float X_diff = X_round - X_cen;
//     float Y_diff = Y_round - Y_cen;
//     float Z_diff = Z_round - Z_cen;

//     int type_k = atoms[ii]-1;  // ensure consistent encoding
//     float h_k  = h[type_k];
//     float b_k  = b[type_k];

//     float sum_z_exp = 0.0f, sum_z2_exp=0.0f;
//     float sum_Dz_exp_l2_z_b = 0.0f;
//     for (int w = -halfWidth; w <= halfWidth; ++w) {
//         float dz         = static_cast<float>(w) + Z_diff;
//         float exp_l2_z_b = __expf(-dz*dz * b_k);
//         sum_z_exp     += exp_l2_z_b;        //__expf, expf
//         sum_z2_exp    += dz*dz*exp_l2_z_b;
//         sum_Dz_exp_l2_z_b += dz*exp_l2_z_b;
//         // if(jj==1 && ii == 9841){
//         //     printf("%.3f,", exp_l2_z_b);
//         //     if(w==halfWidth) printf("\n");
//         // }        
//     }
//     // (1,9841) (4,1836)
//     // if(jj==1 && ii == 9841){
//     //     // printf("%.4f\n", b_k);
//     //     printf("%.4f, %.4f\n", X_round, X_cen);
//     // }

//     for (int u = -halfWidth; u <= halfWidth; ++u) {
//         for (int v = -halfWidth; v <= halfWidth; ++v) {
//             float dx = static_cast<float>(u) + X_diff;
//             float dy = static_cast<float>(v) + Y_diff;
//             float l2_xy = dx*dx + dy*dy;
//             float exp_l2_xy_b = __expf(-l2_xy * b_k);

//             float p_ii   = h_k * exp_l2_xy_b * sum_z_exp;
//             float p_ii_b = p_ii*b_k; 
            

//             float grad_exp = exp_l2_xy_b * sum_z2_exp;
//             float b_ii     = h_k*grad_exp + p_ii*l2_xy;

//             float R2_Dx = ( Rs_jj[0]*dx + Rs_jj[1]*dy ) * p_ii_b;
//             float R2_Dy = ( Rs_jj[3]*dx + Rs_jj[4]*dy ) * p_ii_b;
//             float R2_Dz = ( Rs_jj[6]*dx + Rs_jj[7]*dy ) * p_ii_b;

//             float sum_Dz_exp = sum_Dz_exp_l2_z_b * exp_l2_xy_b;
//             float sum_Dz_hb  = sum_Dz_exp * h_k * b_k;

//             float x_ii = R2_Dx + Rs_jj[2]*sum_Dz_hb;
//             float y_ii = R2_Dy + Rs_jj[5]*sum_Dz_hb;
//             float z_ii = R2_Dz + Rs_jj[8]*sum_Dz_hb;
            
//             // if(jj==1 && ii == 9841){
//             //     printf("%.4f, ", dx);
//             //     if(v==halfWidth) printf("\n");
//             // }
//             // int indx = X_round + u + (N1/2);
//             // int indy = Y_round + v + (N2/2);
//             int indx = X_round_int + u + N1_half;
//             int indy = Y_round_int + v + N2_half;
//             if (indx >= 0 && indx < N1 && indy >= 0 && indy < N2) {
//                 // use 64-bit for index math
//                 size_t proj_idx = static_cast<size_t>(jj)*N1*N2 + indy*N1 + indx;
//                 atomicAdd(&proj[proj_idx], p_ii);

//                 size_t grad_idx = static_cast<size_t>(ii)*num_pj*b_size*b_size + jj*b_size*b_size + (v+halfWidth)*b_size + u+halfWidth;
//                 // Grad_h[u,v,ii,jj] = p_ii;
//                 // Grad_b[u,v,ii,jj] = b_ii;
//                 // Grad_x[u,v,ii,jj] = x_ii;
//                 // Grad_y[u,v,ii,jj] = y_ii;
//                 // Grad_z[u,v,ii,jj] = z_ii;
//                 // Grad_h[grad_idx] = p_ii;
//                 // Grad_b[grad_idx] = b_ii;
//                 Grad_x[grad_idx] = x_ii;
//                 Grad_y[grad_idx] = y_ii;
//                 Grad_z[grad_idx] = z_ii;                
//             }
//         }
//     }
// }
//     // // Accumulate into Grad
//     // for (int u = -halfWidth; u <= halfWidth; ++u) {
//     //     for (int v = -halfWidth; v <= halfWidth; ++v) {
//     //         int indx = X_round + u + (N1/2);
//     //         int indy = Y_round + v + (N2/2);
//     //         int grad_idx = pj*N1*N2 + N1*indy + indx;
//     //         // grad:[N1,N2,num_pj]
//     //         atomicAdd(&Grad[grad_idx], pj_j_h[u + halfWidth][v + halfWidth]);
//     //     }
//     // }



// // extern __global__ void accumulateProjectionKernel(
// //     const float* X_rot, const float* Y_rot, const float* Z_rot,
// //     const int* atom, const int* num_atom_type,
// //     const float* X_crop, const float* Y_crop, const float* Z_crop,
// //     const float* b, const float* h,
// //     float* Grad,
// //     int N1, int N2, int halfWidth,
// //     int num_pj, int atom_type_num,
// //     int max_atoms_per_type);


// __global__ void backward_kernel(
//     float* model,
//     const float* Rs,
//     const int* atoms, int num_atoms,
//     const float* b, const float* h,
//     const float* res,          // <-- use proj instead of undeclared res
//     // float* Grad_h,
//     // float* Grad_b,    
//     const float* Grad_x,
//     const float* Grad_y,
//     const float* Grad_z,
//     float* Grad_model,
//     float dt,
//     int N1, int N2, int halfWidth,
//     int num_pj)
// {
//     int jj = blockIdx.x;  
//     int ii = blockIdx.y * blockDim.y + threadIdx.y;  

//     if (jj >= num_pj || ii >= num_atoms) return;

//     int b_size = 2*halfWidth+1;
//     const float* Rs_jj = Rs + 9*jj;
//     float step_sz = 1.0f;
//     // float dt = step_sz / (halfWidth*halfWidth*num_pj*N1*N2*h[0]*b[0]);

//     float X_cen = Rs_jj[0]*model[3*ii] + Rs_jj[3]*model[3*ii+1] + Rs_jj[6]*model[3*ii+2];
//     float Y_cen = Rs_jj[1]*model[3*ii] + Rs_jj[4]*model[3*ii+1] + Rs_jj[7]*model[3*ii+2];
//     float Z_cen = Rs_jj[2]*model[3*ii] + Rs_jj[5]*model[3*ii+1] + Rs_jj[8]*model[3*ii+2];

//     // int X_round_int = __float2int_rn(X_cen);
//     // int Y_round_int = __float2int_rn(Y_cen);
//     // int Z_round_int = __float2int_rn(Z_cen);
//     float X_round = roundf(X_cen);
//     float Y_round = roundf(Y_cen);
//     float Z_round = roundf(Z_cen);    
//     int X_round_int = static_cast<int>(X_round);
//     int Y_round_int = static_cast<int>(Y_round);

//     int N1_half = N1/2;
//     int N2_half = N2/2;

//     float grad_ii_x = 0.0f;
//     float grad_ii_y = 0.0f;
//     float grad_ii_z = 0.0f;

//     for (int u = -halfWidth; u <= halfWidth; ++u) {
//         for (int v = -halfWidth; v <= halfWidth; ++v) {
//             int indx = X_round_int + u + N1_half;
//             int indy = Y_round_int + v + N2_half;

//             if (indx < 0 || indx >= N1 || indy < 0 || indy >= N2) continue;
            
//             size_t res_idx = ((size_t)jj)*N1*N2 + (size_t)indy*N1 + (size_t)indx;
//             size_t grad_idx = (size_t)ii*num_pj*b_size*b_size
//                             + jj*b_size*b_size
//                             + (v+halfWidth)*b_size
//                             + (u+halfWidth);

//             grad_ii_x += res[res_idx] * Grad_x[grad_idx];
//             grad_ii_y += res[res_idx] * Grad_y[grad_idx];
//             grad_ii_z += res[res_idx] * Grad_z[grad_idx];
//         }
//     }
//     atomicAdd(&Grad_model[3*ii]  , dt*grad_ii_x);
//     atomicAdd(&Grad_model[3*ii+1], dt*grad_ii_y);
//     atomicAdd(&Grad_model[3*ii+2], dt*grad_ii_z);
//     // atomicAdd(&model[3*ii]  , -dt*grad_ii_x);
//     // atomicAdd(&model[3*ii+1], -dt*grad_ii_y);
//     // atomicAdd(&model[3*ii+2], -dt*grad_ii_z);    
// }


// void __global__ diff_kernel(
//     const float* AA,
//     const float* BB,
//     float* diff,
//     size_t N){
//     const size_t ii = blockDim.x * blockIdx.x + threadIdx.x;

//     if (ii<N) {
//         diff[ii] = AA[ii] - BB[ii];
//     }
// }   

// void __global__ update_model_kernel(
//     float* model,
//     const float* model_in,
//     float* Grad_model,
//     const float scale,
//     size_t num_atoms){
//     const size_t ii = blockDim.x * blockIdx.x + threadIdx.x;

//     // 0. out of range
//     if (ii>=num_atoms) return;
    
//     // 1. ii<num_atoms
//     // float grad_x =   - model[3*ii]   + model_in[3*ii];
//     // float grad_y =   - model[3*ii+1] + model_in[3*ii+1];
//     // float grad_z =   - model[3*ii+2] + model_in[3*ii+2];
//     float grad_x =   Grad_model[3*ii]   - model[3*ii]   + model_in[3*ii];
//     float grad_y =   Grad_model[3*ii+1] - model[3*ii+1] + model_in[3*ii+1];
//     float grad_z =   Grad_model[3*ii+2] - model[3*ii+2] + model_in[3*ii+2];
//     float diff_norm = sqrtf( grad_x*grad_x + grad_y*grad_y + grad_z*grad_z );

//     if(diff_norm>scale){
//         grad_x = grad_x*scale/diff_norm;
//         grad_y = grad_y*scale/diff_norm;
//         grad_z = grad_z*scale/diff_norm;
//         model[3*ii]   = model_in[3*ii]   - grad_x;
//         model[3*ii+1] = model_in[3*ii+1] - grad_y;
//         model[3*ii+2] = model_in[3*ii+2] - grad_z;
//     }
//     else{
//         model[3*ii]   -= grad_x;
//         model[3*ii+1] -= grad_y;
//         model[3*ii+2] -= grad_z;
//     }
// }   


// forward_kernel, diff_kernel, backward_kernel, normalized_grad_kernel must be declared elsewhere
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    mxInitGPU();
    const int num_threads = 256;

    // === [1] Parse inputs ===
    const mxGPUArray* model_in_gpu = mxGPUCreateFromMxArray(prhs[0]); // [3, num_atoms]
    const mxGPUArray* atom_gpu  = mxGPUCreateFromMxArray(prhs[1]); // [num_atoms]
    const mxGPUArray* mea_proj_gpu = mxGPUCreateFromMxArray(prhs[2]); // measured projection [N1,N2,num_pj]
    const mxGPUArray* fixedfa_gpu  = mxGPUCreateFromMxArray(prhs[3]); // measured projection [N1,N2,num_pj]
    // const mxArray* fixedfa = prhs[3];
    const mxGPUArray* b_gpu     = mxGPUCreateFromMxArray(prhs[4]); // [num_types]
    const mxGPUArray* h_gpu     = mxGPUCreateFromMxArray(prhs[5]); // [num_types]
    const mxGPUArray* Rs_gpu    = mxGPUCreateFromMxArray(prhs[6]); // [3,3,num_pj]

    float scale    = static_cast<float>(mxGetScalar(prhs[7]));
    int iterations = static_cast<int>(mxGetScalar(prhs[8]));
    float dt       = static_cast<float>(mxGetScalar(prhs[9]));
    int halfWidth  = static_cast<int>(mxGetScalar(prhs[10]));
    int b_size     = 2*halfWidth+1;
    std::cout << "iterations = " << iterations <<std::endl;

    // === [2] Dimension checks ===
    // check input: model [3,num_atoms]
    const mwSize* modelDims = mxGPUGetDimensions(model_in_gpu);
    if (mxGPUGetNumberOfDimensions(model_in_gpu) != 2 || modelDims[0] != 3) {
        mexErrMsgIdAndTxt("invalidModel","model must be [3,num_atoms]");
    }
    if (mxGPUGetClassID(model_in_gpu) != mxSINGLE_CLASS)
        mexErrMsgIdAndTxt("typeError","model must be single precision");
    int num_atoms = static_cast<int>(modelDims[1]);


    // check input: atoms type [3,num_atoms]
    int atoms_ndims = mxGPUGetNumberOfDimensions(atom_gpu);
    const mwSize* atomDims = mxGPUGetDimensions(atom_gpu);
    if (atoms_ndims != 2 ||
        !((atomDims[0] == 1 && atomDims[1] == num_atoms) ||
          (atomDims[0] == num_atoms && atomDims[1] == 1))) {
        mexErrMsgIdAndTxt("invalidAtoms","atoms must be [1,num_atoms] or [num_atoms,1]");
    }
    if (mxGPUGetClassID(atom_gpu) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("typeError","atoms must be int32");


    // check input: projeciton [N1,N2,num_pj]
    const mwSize* pjDims = mxGPUGetDimensions(mea_proj_gpu);
    if (mxGPUGetNumberOfDimensions(mea_proj_gpu) != 3) {
        mexErrMsgIdAndTxt("invalidProj","projection must be [N1,N2,num_pj]");
    }
    if (mxGPUGetClassID(mea_proj_gpu) != mxSINGLE_CLASS)
        mexErrMsgIdAndTxt("typeError","projection must be single precision");
    const int N1 = pjDims[0];
    const int N2 = pjDims[1];
    const int num_pj = pjDims[2];

    // check input: h, b
    const mwSize* bDims = mxGPUGetDimensions(b_gpu);
    const mwSize* hDims = mxGPUGetDimensions(h_gpu);
    mwSize bNumel = 1; for (int d=0; d<mxGPUGetNumberOfDimensions(b_gpu); ++d) bNumel *= bDims[d];
    mwSize hNumel = 1; for (int d=0; d<mxGPUGetNumberOfDimensions(h_gpu); ++d) hNumel *= hDims[d];
    if (bNumel != hNumel) {
        mexErrMsgIdAndTxt("invalidBH","b and h must be vectors of equal length");
    }
    int num_types = static_cast<int>(bNumel);

    // check input: Rs [3,3,num_pj]
    const mwSize* RsDims = mxGPUGetDimensions(Rs_gpu);
    if (mxGPUGetNumberOfDimensions(Rs_gpu)!=3 || RsDims[0]!=3 || RsDims[1]!=3 || RsDims[2]!=num_pj) {
        mexErrMsgIdAndTxt("invalidRs","Rs must be [3,3,num_pj]");
    }
    if (mxGPUGetClassID(Rs_gpu) != mxSINGLE_CLASS)
        mexErrMsgIdAndTxt("typeError","Rs must be single precision");    
    
    size_t num_pts = static_cast<size_t>(N1) * N2 * num_pj;

    // // === [3] Get device pointers ===
    // const float* model_in = static_cast<const float*>(mxGPUGetDataReadOnly(model_in_gpu));
    // const int*   atoms    = static_cast<const int*>(mxGPUGetDataReadOnly(atom_gpu));
    // const float* mea_proj = static_cast<const float*>(mxGPUGetDataReadOnly(mea_proj_gpu));
    // const float* b        = static_cast<const float*>(mxGPUGetDataReadOnly(b_gpu));
    // const float* h        = static_cast<const float*>(mxGPUGetDataReadOnly(h_gpu));
    // const float* Rs       = static_cast<const float*>(mxGPUGetDataReadOnly(Rs_gpu));

    // // Allocate a new GPU array for a writable copy
    // // mxGPUArray* model_copy_gpu = mxGPUCreateGPUArray(
    // //     mxGPUGetNumberOfDimensions(model_in_gpu),
    // //     modelDims, mxSINGLE_CLASS, mxREAL, MX_GPU_DO_NOT_INITIALIZE
    // // );
    // // float* model = static_cast<float*>(mxGPUGetData(model_copy_gpu));
    // float* model;
    // cudaMalloc( &model, 3*num_atoms*sizeof(float) );        
    // cudaMemcpy(model, model_in, 3*num_atoms * sizeof(float), cudaMemcpyDeviceToDevice);

    // // === [4] Allocate output projection ===
    // // mxGPUArray* Grad_gpu = mxGPUCreateGPUArray(3, pjDims, mxSINGLE_CLASS, mxREAL, MX_GPU_DO_NOT_INITIALIZE);
    // // float* proj = static_cast<float*>(mxGPUGetData(Grad_gpu));
    // float * proj;
    // cudaMalloc( &proj, num_pts*sizeof(float) );    
    // // cudaMemset(proj, 0, num_pts * sizeof(float));

    // // === [5] Launch kernels ===
    // // a. forward projection
    // dim3 grid(num_pj,(num_atoms + num_threads - 1) / num_threads);
    // dim3 block(1, num_threads);
    // dim3 block0(num_threads);

    size_t num_b_pts = (size_t)b_size*b_size*num_atoms*num_pj;
    mwSize grad_x_Dims[4] = {(size_t)b_size, (size_t)b_size, (size_t)num_pj, (size_t)num_atoms, };
    // // mxGPUArray* Grad_x_gpu = mxGPUCreateGPUArray(4, grad_x_Dims, mxSINGLE_CLASS, mxREAL, MX_GPU_DO_NOT_INITIALIZE);
    // // float* Grad_x = static_cast<float*>(mxGPUGetData(Grad_x_gpu));
    // float * Grad_x;
    // cudaMalloc( &Grad_x, num_b_pts*sizeof(float) );       
    // // cudaMemset(Grad_x, 0, num_b_pts * sizeof(float));

    // // mxGPUArray* Grad_y_gpu = mxGPUCreateGPUArray(4, grad_x_Dims, mxSINGLE_CLASS, mxREAL, MX_GPU_DO_NOT_INITIALIZE);
    // // float* Grad_y = static_cast<float*>(mxGPUGetData(Grad_y_gpu));
    // float * Grad_y;
    // cudaMalloc( &Grad_y, num_b_pts*sizeof(float) );           
    // // cudaMemset(Grad_y, 0, num_b_pts * sizeof(float));

    // // mxGPUArray* Grad_z_gpu = mxGPUCreateGPUArray(4, grad_x_Dims, mxSINGLE_CLASS, mxREAL, MX_GPU_DO_NOT_INITIALIZE);
    // // float* Grad_z = static_cast<float*>(mxGPUGetData(Grad_z_gpu));
    // float * Grad_z;
    // cudaMalloc( &Grad_z, num_b_pts*sizeof(float) );           
    // // cudaMemset(Grad_z, 0, num_b_pts * sizeof(float));

    // dim3 grid_diff((num_pts + num_threads - 1) / num_threads);

    // float * proj_diff;
    // cudaMalloc( &proj_diff, num_pts*sizeof(float) );           
    // // cudaMemset(proj_diff, 0, num_pts * sizeof(float));
    
    // // FFTFilterR2C filter(fixedfa, num_pj);
    // FFTFilterR2C filter(fixedfa_gpu, num_pj);

    // float * Grad_model;
    // cudaMalloc( &Grad_model, 3*num_atoms*sizeof(float) );      
    // // cudaMemset(Grad_model, 0, 3*num_atoms * sizeof(float));
        
    // for(int iter=0; iter<iterations; iter++){
    //     // 1: forward projection
    //     cudaMemset(Grad_x, 0, num_b_pts * sizeof(float));
    //     cudaMemset(Grad_y, 0, num_b_pts * sizeof(float));
    //     cudaMemset(Grad_z, 0, num_b_pts * sizeof(float));
    //     cudaMemset(proj, 0, num_pts * sizeof(float));
    //     forward_kernel<<<grid, block>>>(model, Rs, atoms, num_atoms, b, h,
    //                     proj, Grad_x, Grad_y, Grad_z,
    //                     N1, N2, halfWidth, num_pj);
    //     CUDA_CHECK(cudaGetLastError());

    //     // 2. projection filter
    //     filter.apply(proj, proj);

    //     // mxGPUArray* diff_gpu = mxGPUCreateGPUArray(3, pjDims, mxSINGLE_CLASS, mxREAL, MX_GPU_DO_NOT_INITIALIZE);
    //     // float* proj_diff = static_cast<float*>(mxGPUGetData(diff_gpu));

    //     // 3. difference       
    //     // cudaMemset(proj_diff, 0, num_pts * sizeof(float)); 
    //     diff_kernel<<<grid_diff, block0>>>(proj, mea_proj, proj_diff, num_pts);
    //     CUDA_CHECK(cudaGetLastError());

    //     // mxGPUArray* Grad_model_gpu = mxGPUCreateGPUArray(2, modelDims, mxSINGLE_CLASS, mxREAL, MX_GPU_DO_NOT_INITIALIZE);
    //     // float* Grad_model = static_cast<float*>(mxGPUGetData(Grad_model_gpu));

    //     // c. backward
    //     cudaMemset(Grad_model, 0, 3*num_atoms * sizeof(float));
    //     backward_kernel<<<grid, block>>>(model, Rs, atoms, num_atoms, b, h,
    //                                     proj_diff, Grad_x, Grad_y, Grad_z,
    //                                     Grad_model, dt, N1, N2, halfWidth, num_pj);
    //     CUDA_CHECK(cudaGetLastError());


    //     // d. normalization
    //     // float scale = 1.0f;
    //     dim3 grid_norm((num_atoms + num_threads - 1) / num_threads);
    //     update_model_kernel<<<grid_norm, block0>>>(model, model_in, Grad_model, scale, num_atoms);
    //     CUDA_CHECK(cudaGetLastError());

    //     cudaDeviceSynchronize();
    // }

    TomoModelOptimizer optimizer(
        model_in_gpu, atom_gpu, mea_proj_gpu, b_gpu, h_gpu, Rs_gpu, fixedfa_gpu,
        num_atoms, num_pj, num_pts, b_size,
         N1,  N2, halfWidth, 
        num_threads, dt, scale
    );
    optimizer.run(iterations);
    float* model    = optimizer.get_model();
    float* Grad_x   = optimizer.get_Grad_x();
    float* Grad_y   = optimizer.get_Grad_y();
    float* Grad_z   = optimizer.get_Grad_z();
    float* Grad_model   = optimizer.get_Grad_model();
    float* proj   = optimizer.get_proj();
    float* proj_diff   = optimizer.get_proj_diff();

    // === [6] Return output ===
    // plhs[0] = mxGPUCreateMxArrayOnGPU(model_copy_gpu);
    // plhs[1] = mxGPUCreateMxArrayOnGPU(Grad_gpu);
    // plhs[2] = mxGPUCreateMxArrayOnGPU(Grad_x_gpu);
    // plhs[3] = mxGPUCreateMxArrayOnGPU(Grad_y_gpu);
    // plhs[4] = mxGPUCreateMxArrayOnGPU(Grad_z_gpu);
    // plhs[5] = mxGPUCreateMxArrayOnGPU(Grad_model_gpu);

    // mxGPUDestroyGPUArray(model_copy_gpu);
    // mxGPUDestroyGPUArray(Grad_gpu);
    // mxGPUDestroyGPUArray(Grad_x_gpu);
    // mxGPUDestroyGPUArray(Grad_y_gpu);
    // mxGPUDestroyGPUArray(Grad_z_gpu);
    // mxGPUDestroyGPUArray(Grad_model_gpu);
    // mxGPUDestroyGPUArray(diff_gpu);

    // === [6] Copy results back to MATLAB ===
    // model copy
    mwSize outDimsModel[2] = {3, (mwSize)num_atoms};
    plhs[0] = mxCreateNumericArray(2, outDimsModel, mxSINGLE_CLASS, mxREAL);
    float* model_cpu = (float*)mxGetData(plhs[0]);
    cudaMemcpy(model_cpu, model, 3*num_atoms*sizeof(float), cudaMemcpyDeviceToHost);

    // proj
    plhs[1] = mxCreateNumericArray(3, pjDims, mxSINGLE_CLASS, mxREAL);
    float* proj_cpu = (float*)mxGetData(plhs[1]);
    cudaMemcpy(proj_cpu, proj, num_pts*sizeof(float), cudaMemcpyDeviceToHost);

    // Grad_x
    plhs[2] = mxCreateNumericArray(4, grad_x_Dims, mxSINGLE_CLASS, mxREAL);
    float* Grad_x_cpu = (float*)mxGetData(plhs[2]);
    cudaMemcpy(Grad_x_cpu, Grad_x, num_b_pts*sizeof(float), cudaMemcpyDeviceToHost);

    // Grad_y
    plhs[3] = mxCreateNumericArray(4, grad_x_Dims, mxSINGLE_CLASS, mxREAL);
    float* Grad_y_cpu = (float*)mxGetData(plhs[3]);
    cudaMemcpy(Grad_y_cpu, Grad_y, num_b_pts*sizeof(float), cudaMemcpyDeviceToHost);

    // Grad_z
    plhs[4] = mxCreateNumericArray(4, grad_x_Dims, mxSINGLE_CLASS, mxREAL);
    float* Grad_z_cpu = (float*)mxGetData(plhs[4]);
    cudaMemcpy(Grad_z_cpu, Grad_z, num_b_pts*sizeof(float), cudaMemcpyDeviceToHost);

    // Grad_model
    // plhs[5] = mxCreateNumericArray(2, modelDims, mxSINGLE_CLASS, mxREAL);
    // float* Grad_model_cpu = (float*)mxGetData(plhs[5]);
    // cudaMemcpy(Grad_model_cpu, Grad_model, 3*num_atoms*sizeof(float), cudaMemcpyDeviceToHost);


    // === [7] Cleanup (do not destroy arrays returned to MATLAB) ===
    mxGPUDestroyGPUArray(model_in_gpu);
    mxGPUDestroyGPUArray(atom_gpu);
    mxGPUDestroyGPUArray(mea_proj_gpu);
    mxGPUDestroyGPUArray(b_gpu);
    mxGPUDestroyGPUArray(h_gpu);
    mxGPUDestroyGPUArray(Rs_gpu);
    mxGPUDestroyGPUArray(fixedfa_gpu);

    // cudaFree(model);
    // cudaFree(proj);
    // cudaFree(Grad_x);
    // cudaFree(Grad_y);
    // cudaFree(Grad_z);
    // cudaFree(Grad_model);
    // cudaFree(proj_diff);

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        mexErrMsgIdAndTxt("accumulateProjection:cudaError", cudaGetErrorString(err));
    }    
}

