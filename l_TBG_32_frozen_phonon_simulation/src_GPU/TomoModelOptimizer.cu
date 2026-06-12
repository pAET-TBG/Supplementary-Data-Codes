#include "mex.h"
#include "omp.h"
#include <cmath>
#include "gpu/mxGPUArray.h"
#include <cuda_runtime.h>
#include <cufft.h>
#include "FFTFilterR2C.cu"

#ifndef CUDA_CHECK_DEFINED
#define CUDA_CHECK_DEFINED
#define CUDA_CHECK(call) do { \
    cudaError_t err = (call); \
    if (err != cudaSuccess) { \
        mexErrMsgIdAndTxt("cudaError", cudaGetErrorString(err)); \
    } \
} while(0)
#endif // CUDA_CHECK_DEFINED

__global__ void forward_kernel(
    const float* model,
    const float* Rs,
    const int* atoms, int num_atoms,
    const float* b, const float* h,
    float* proj,
    // float* Grad_h,
    // float* Grad_b,
    float* Grad_x,
    float* Grad_y,
    float* Grad_z,
    int N1, int N2, int halfWidth,
    int num_pj)
{
    int jj = blockIdx.x;  
    int ii = blockIdx.y * blockDim.y + threadIdx.y;  

    if (jj >= num_pj || ii >= num_atoms) return;


    const float* Rs_jj = Rs + 9*jj;

    float X_cen = Rs_jj[0]*model[3*ii] + Rs_jj[3]*model[3*ii+1] + Rs_jj[6]*model[3*ii+2];
    float Y_cen = Rs_jj[1]*model[3*ii] + Rs_jj[4]*model[3*ii+1] + Rs_jj[7]*model[3*ii+2];
    float Z_cen = Rs_jj[2]*model[3*ii] + Rs_jj[5]*model[3*ii+1] + Rs_jj[8]*model[3*ii+2];

    // int X_round = __float2int_rn(X_cen);
    // int Y_round = __float2int_rn(Y_cen);
    // int Z_round = __float2int_rn(Z_cen);

    float X_round = roundf(X_cen);
    float Y_round = roundf(Y_cen);
    float Z_round = roundf(Z_cen);    

    // float X_round = roundAwayFromZero(X_cen);
    // float Y_round = roundAwayFromZero(Y_cen);
    // float Z_round = roundAwayFromZero(Z_cen);     
    
    int X_round_int = X_round;
    int Y_round_int = Y_round;
    // int Z_round_int = static_cast<int>(Z_round);

    // if(jj==1 && ii == 9841){
    //     printf("%d,%d: %f,%f,%f\n",jj,ii, X_cen, Y_cen, Z_cen);
    //     printf("%d,%d: %f,%f,%f\n",jj,ii, X_round, Y_round, Z_round);
    //     printf("%d,%d: %d,%d\n",jj,ii, X_round_int, Y_round_int);
    // }

    int N1_half = N1/2;
    int N2_half = N2/2;

    int b_size = 2*halfWidth+1;

    float X_diff = X_round - X_cen;
    float Y_diff = Y_round - Y_cen;
    float Z_diff = Z_round - Z_cen;

    int type_k = atoms[ii]-1;  // ensure consistent encoding
    float h_k  = h[type_k];
    float b_k  = b[type_k];

    float sum_z_exp = 0.0f, sum_z2_exp=0.0f;
    float sum_Dz_exp_l2_z_b = 0.0f;
    for (int w = -halfWidth; w <= halfWidth; ++w) {
        float dz         = static_cast<float>(w) + Z_diff;
        float exp_l2_z_b = __expf(-dz*dz * b_k);
        sum_z_exp     += exp_l2_z_b;        //__expf, expf
        sum_z2_exp    += dz*dz*exp_l2_z_b;
        sum_Dz_exp_l2_z_b += dz*exp_l2_z_b;
        // if(jj==1 && ii == 9841){
        //     printf("%.3f,", exp_l2_z_b);
        //     if(w==halfWidth) printf("\n");
        // }        
    }
    // (1,9841) (4,1836)
    // if(jj==1 && ii == 9841){
    //     // printf("%.4f\n", b_k);
    //     printf("%.4f, %.4f\n", X_round, X_cen);
    // }

    for (int u = -halfWidth; u <= halfWidth; ++u) {
        for (int v = -halfWidth; v <= halfWidth; ++v) {
            float dx = static_cast<float>(u) + X_diff;
            float dy = static_cast<float>(v) + Y_diff;
            float l2_xy = dx*dx + dy*dy;
            float exp_l2_xy_b = __expf(-l2_xy * b_k);

            float p_ii   = h_k * exp_l2_xy_b * sum_z_exp;
            float p_ii_b = p_ii*b_k; 
            

            float grad_exp = exp_l2_xy_b * sum_z2_exp;
            float b_ii     = h_k*grad_exp + p_ii*l2_xy;

            float R2_Dx = ( Rs_jj[0]*dx + Rs_jj[1]*dy ) * p_ii_b;
            float R2_Dy = ( Rs_jj[3]*dx + Rs_jj[4]*dy ) * p_ii_b;
            float R2_Dz = ( Rs_jj[6]*dx + Rs_jj[7]*dy ) * p_ii_b;

            float sum_Dz_exp = sum_Dz_exp_l2_z_b * exp_l2_xy_b;
            float sum_Dz_hb  = sum_Dz_exp * h_k * b_k;

            float x_ii = R2_Dx + Rs_jj[2]*sum_Dz_hb;
            float y_ii = R2_Dy + Rs_jj[5]*sum_Dz_hb;
            float z_ii = R2_Dz + Rs_jj[8]*sum_Dz_hb;
            
            // if(jj==1 && ii == 9841){
            //     printf("%.4f, ", dx);
            //     if(v==halfWidth) printf("\n");
            // }
            // int indx = X_round + u + (N1/2);
            // int indy = Y_round + v + (N2/2);
            int indx = X_round_int + u + N1_half;
            int indy = Y_round_int + v + N2_half;
            if (indx >= 0 && indx < N1 && indy >= 0 && indy < N2) {
                // use 64-bit for index math
                size_t proj_idx = static_cast<size_t>(jj)*N1*N2 + indy*N1 + indx;
                atomicAdd(&proj[proj_idx], p_ii);

                size_t grad_idx = static_cast<size_t>(ii)*num_pj*b_size*b_size + jj*b_size*b_size + (v+halfWidth)*b_size + u+halfWidth;
                // Grad_h[u,v,ii,jj] = p_ii;
                // Grad_b[u,v,ii,jj] = b_ii;
                // Grad_x[u,v,ii,jj] = x_ii;
                // Grad_y[u,v,ii,jj] = y_ii;
                // Grad_z[u,v,ii,jj] = z_ii;
                // Grad_h[grad_idx] = p_ii;
                // Grad_b[grad_idx] = b_ii;
                Grad_x[grad_idx] = x_ii;
                Grad_y[grad_idx] = y_ii;
                Grad_z[grad_idx] = z_ii;                
            }
        }
    }
}
    // // Accumulate into Grad
    // for (int u = -halfWidth; u <= halfWidth; ++u) {
    //     for (int v = -halfWidth; v <= halfWidth; ++v) {
    //         int indx = X_round + u + (N1/2);
    //         int indy = Y_round + v + (N2/2);
    //         int grad_idx = pj*N1*N2 + N1*indy + indx;
    //         // grad:[N1,N2,num_pj]
    //         atomicAdd(&Grad[grad_idx], pj_j_h[u + halfWidth][v + halfWidth]);
    //     }
    // }



// extern __global__ void accumulateProjectionKernel(
//     const float* X_rot, const float* Y_rot, const float* Z_rot,
//     const int* atom, const int* num_atom_type,
//     const float* X_crop, const float* Y_crop, const float* Z_crop,
//     const float* b, const float* h,
//     float* Grad,
//     int N1, int N2, int halfWidth,
//     int num_pj, int atom_type_num,
//     int max_atoms_per_type);


__global__ void backward_kernel(
    float* model,
    const float* Rs,
    const int* atoms, int num_atoms,
    const float* b, const float* h,
    const float* res,          // <-- use proj instead of undeclared res
    // float* Grad_h,
    // float* Grad_b,    
    const float* Grad_x,
    const float* Grad_y,
    const float* Grad_z,
    float* Grad_model,
    float dt,
    int N1, int N2, int halfWidth,
    int num_pj)
{
    int jj = blockIdx.x;  
    int ii = blockIdx.y * blockDim.y + threadIdx.y;  

    if (jj >= num_pj || ii >= num_atoms) return;

    int b_size = 2*halfWidth+1;
    const float* Rs_jj = Rs + 9*jj;
    float step_sz = 1.0f;
    // float dt = step_sz / (halfWidth*halfWidth*num_pj*N1*N2*h[0]*b[0]);

    float X_cen = Rs_jj[0]*model[3*ii] + Rs_jj[3]*model[3*ii+1] + Rs_jj[6]*model[3*ii+2];
    float Y_cen = Rs_jj[1]*model[3*ii] + Rs_jj[4]*model[3*ii+1] + Rs_jj[7]*model[3*ii+2];
    float Z_cen = Rs_jj[2]*model[3*ii] + Rs_jj[5]*model[3*ii+1] + Rs_jj[8]*model[3*ii+2];

    // int X_round_int = __float2int_rn(X_cen);
    // int Y_round_int = __float2int_rn(Y_cen);
    // int Z_round_int = __float2int_rn(Z_cen);
    float X_round = roundf(X_cen);
    float Y_round = roundf(Y_cen);
    float Z_round = roundf(Z_cen);    
    int X_round_int = static_cast<int>(X_round);
    int Y_round_int = static_cast<int>(Y_round);

    int N1_half = N1/2;
    int N2_half = N2/2;

    float grad_ii_x = 0.0f;
    float grad_ii_y = 0.0f;
    float grad_ii_z = 0.0f;

    for (int u = -halfWidth; u <= halfWidth; ++u) {
        for (int v = -halfWidth; v <= halfWidth; ++v) {
            int indx = X_round_int + u + N1_half;
            int indy = Y_round_int + v + N2_half;

            if (indx < 0 || indx >= N1 || indy < 0 || indy >= N2) continue;
            
            size_t res_idx = ((size_t)jj)*N1*N2 + (size_t)indy*N1 + (size_t)indx;
            size_t grad_idx = (size_t)ii*num_pj*b_size*b_size
                            + jj*b_size*b_size
                            + (v+halfWidth)*b_size
                            + (u+halfWidth);

            grad_ii_x += res[res_idx] * Grad_x[grad_idx];
            grad_ii_y += res[res_idx] * Grad_y[grad_idx];
            grad_ii_z += res[res_idx] * Grad_z[grad_idx];
        }
    }
    atomicAdd(&Grad_model[3*ii]  , dt*grad_ii_x);
    atomicAdd(&Grad_model[3*ii+1], dt*grad_ii_y);
    atomicAdd(&Grad_model[3*ii+2], dt*grad_ii_z);
    // atomicAdd(&model[3*ii]  , -dt*grad_ii_x);
    // atomicAdd(&model[3*ii+1], -dt*grad_ii_y);
    // atomicAdd(&model[3*ii+2], -dt*grad_ii_z);    
}


void __global__ diff_kernel(
    const float* AA,
    const float* BB,
    float* diff,
    size_t N){
    const size_t ii = blockDim.x * blockIdx.x + threadIdx.x;

    if (ii<N) {
        diff[ii] = AA[ii] - BB[ii];
    }
}   

void __global__ update_model_kernel(
    float* model,
    const float* model_in,
    float* Grad_model,
    const float scale,
    size_t num_atoms){
    const size_t ii = blockDim.x * blockIdx.x + threadIdx.x;

    // 0. out of range
    if (ii>=num_atoms) return;
    
    // 1. ii<num_atoms
    // ||Gx - xk + x0||^2 = R^2
    // ||G||^2*x^2 - 2x<G,xk-x0> + ||xk-x0||^2-R^2 = 0
    float aa     = Grad_model[3*ii]*Grad_model[3*ii] + Grad_model[3*ii+1]*Grad_model[3*ii+1] + Grad_model[3*ii+2]*Grad_model[3*ii+2];
    float diff_x = model[3*ii]   - model_in[3*ii];
    float diff_y = model[3*ii+1] - model_in[3*ii+1];
    float diff_z = model[3*ii+2] - model_in[3*ii+2];
    float bb = Grad_model[3*ii]*diff_x + Grad_model[3*ii+1]*diff_y + Grad_model[3*ii+2]*diff_z;
    float cc = diff_x*diff_x + diff_y*diff_y + diff_z*diff_z - scale*scale;
    cc = fminf(cc,0);
    // solve a*x^2 - 2b*x + c = 0;
    float sqrt_delta = sqrtf( bb*bb - aa*cc );
    // float x1 = (bb - sqrt_delta) / aa; 
    // float x2 = (bb + sqrt_delta) / aa; 
    // float soln = fmaxf(x1,x2); 
    float soln = (bb + sqrt_delta) / aa; 
    soln = fminf(soln, 1);
    model[3*ii]   -= soln*Grad_model[3*ii];
    model[3*ii+1] -= soln*Grad_model[3*ii+1];
    model[3*ii+2] -= soln*Grad_model[3*ii+2];

    // float grad_x =   - model[3*ii]   + model_in[3*ii];
    // float grad_y =   - model[3*ii+1] + model_in[3*ii+1];
    // float grad_z =   - model[3*ii+2] + model_in[3*ii+2];
    // float grad_x =   Grad_model[3*ii]   - model[3*ii]   + model_in[3*ii];
    // float grad_y =   Grad_model[3*ii+1] - model[3*ii+1] + model_in[3*ii+1];
    // float grad_z =   Grad_model[3*ii+2] - model[3*ii+2] + model_in[3*ii+2];
    // float diff_norm = sqrtf( grad_x*grad_x + grad_y*grad_y + grad_z*grad_z );

    // if(diff_norm>scale){
    //     grad_x = grad_x*scale/diff_norm;
    //     grad_y = grad_y*scale/diff_norm;
    //     grad_z = grad_z*scale/diff_norm;
    //     model[3*ii]   = model_in[3*ii]   - grad_x;
    //     model[3*ii+1] = model_in[3*ii+1] - grad_y;
    //     model[3*ii+2] = model_in[3*ii+2] - grad_z;
    // }
    // else{
    //     model[3*ii]   -= Grad_model[3*ii];
    //     model[3*ii+1] -= Grad_model[3*ii+1];
    //     model[3*ii+2] -= Grad_model[3*ii+2];
    // }
}   


class TomoModelOptimizer
{
public:
    using Real = float;

    TomoModelOptimizer(
        const mxGPUArray* model_in_gpu,
        const mxGPUArray* atom_gpu,
        const mxGPUArray* mea_proj_gpu,
        const mxGPUArray* b_gpu,
        const mxGPUArray* h_gpu,
        const mxGPUArray* Rs_gpu,
        const mxGPUArray* fixedfa_gpu,
        int num_atoms,
        int num_pj,
        int num_pts,
        int b_size,
        int N1, int N2,
        int halfWidth,
        int num_threads,
        float dt,
        float scale)
        :
        num_atoms_(num_atoms),
        num_pj_(num_pj),
        num_pts_(num_pts),
        b_size_(b_size),
        N1_(N1), N2_(N2),
        halfWidth_(halfWidth),
        num_threads_(num_threads),
        dt_(dt),
        scale_(scale)
    {
        // -----------------------------
        // Store read-only GPU pointers
        // -----------------------------
        model_in_ = static_cast<const Real*>(mxGPUGetDataReadOnly(model_in_gpu));
        atoms_    = static_cast<const int*>(mxGPUGetDataReadOnly(atom_gpu));
        mea_proj_ = static_cast<const Real*>(mxGPUGetDataReadOnly(mea_proj_gpu));
        b_        = static_cast<const Real*>(mxGPUGetDataReadOnly(b_gpu));
        h_        = static_cast<const Real*>(mxGPUGetDataReadOnly(h_gpu));
        Rs_       = static_cast<const Real*>(mxGPUGetDataReadOnly(Rs_gpu));

        // -----------------------------
        // Allocate persistent buffers
        // -----------------------------
        CUDA_CHECK(cudaMalloc(&model_, 3*num_atoms_*sizeof(Real)));
        CUDA_CHECK(cudaMemcpy(model_, model_in_, 3*num_atoms_*sizeof(Real),
                              cudaMemcpyDeviceToDevice));

        CUDA_CHECK(cudaMalloc(&proj_, num_pts_*sizeof(Real)));
        CUDA_CHECK(cudaMalloc(&proj_diff_, num_pts_*sizeof(Real)));

        size_t num_b_pts = (size_t)b_size_*b_size_*num_atoms_*num_pj_;
        CUDA_CHECK(cudaMalloc(&Grad_x_, num_b_pts*sizeof(Real)));
        CUDA_CHECK(cudaMalloc(&Grad_y_, num_b_pts*sizeof(Real)));
        CUDA_CHECK(cudaMalloc(&Grad_z_, num_b_pts*sizeof(Real)));

        CUDA_CHECK(cudaMalloc(&Grad_model_, 3*num_atoms_*sizeof(Real)));

        // -----------------------------
        // Precompute launch geometry
        // -----------------------------
        grid_fwd_  = dim3(num_pj_, (num_atoms_ + num_threads_ - 1) / num_threads_);
        block_fwd_ = dim3(1, num_threads_);
        block1D_   = dim3(num_threads_);
        grid_diff_ = dim3((num_pts_ + num_threads_ - 1) / num_threads_);

        // -----------------------------
        // FFT filter
        // -----------------------------
        // filter_ = std::make_unique<FFTFilterR2C>(fixedfa_gpu, num_pj_);
        filter_ = new FFTFilterR2C(fixedfa_gpu, num_pj_);    
    }

    ~TomoModelOptimizer()
    {
        cudaFree(model_);
        cudaFree(proj_);
        cudaFree(proj_diff_);
        cudaFree(Grad_x_);
        cudaFree(Grad_y_);
        cudaFree(Grad_z_);
        cudaFree(Grad_model_);
        delete filter_;
    }

    // ------------------------------------------------------------
    // Main optimization loop
    // ------------------------------------------------------------
    void run(int iterations)
    {
        for(int iter = 0; iter < iterations; iter++)
        {
            // -----------------------------
            // 1. Forward projection
            // -----------------------------
            cudaMemset(Grad_x_, 0, bBytes());
            cudaMemset(Grad_y_, 0, bBytes());
            cudaMemset(Grad_z_, 0, bBytes());
            cudaMemset(proj_,   0, num_pts_*sizeof(Real));

            forward_kernel<<<grid_fwd_, block_fwd_>>>(
                model_, Rs_, atoms_, num_atoms_, b_, h_,
                proj_, Grad_x_, Grad_y_, Grad_z_,
                N1_, N2_, halfWidth_, num_pj_);
            CUDA_CHECK(cudaGetLastError());

            // -----------------------------
            // 2. FFT filtering
            // -----------------------------
            filter_->apply(proj_, proj_);

            // -----------------------------
            // 3. Difference
            // -----------------------------
            diff_kernel<<<grid_diff_, block1D_>>>(
                proj_, mea_proj_, proj_diff_, num_pts_);
            CUDA_CHECK(cudaGetLastError());

            // -----------------------------
            // 4. Backward projection
            // -----------------------------
            cudaMemset(Grad_model_, 0, 3*num_atoms_*sizeof(Real));

            backward_kernel<<<grid_fwd_, block_fwd_>>>(
                model_, Rs_, atoms_, num_atoms_, b_, h_,
                proj_diff_, Grad_x_, Grad_y_, Grad_z_,
                Grad_model_, dt_,
                N1_, N2_, halfWidth_, num_pj_);
            CUDA_CHECK(cudaGetLastError());

            // -----------------------------
            // 5. Update model
            // -----------------------------
            dim3 grid_norm((num_atoms_ + num_threads_ - 1) / num_threads_);
            update_model_kernel<<<grid_norm, block1D_>>>(
                model_, model_in_, Grad_model_, scale_, num_atoms_);
            CUDA_CHECK(cudaGetLastError());

            cudaDeviceSynchronize();
        }
    }

    float* get_model() const {
        return model_;
    }

    float* get_proj() const {
        return proj_;
    }

    float* get_Grad_x() const {
        return Grad_x_;
    }

    float* get_Grad_y() const {
        return Grad_y_;
    }
    
    float* get_Grad_z() const {
        return Grad_z_;
    }
    
    float* get_Grad_model() const {
        return Grad_model_;
    }    

    float* get_proj_diff() const {
        return proj_diff_;
    }

private:
    // Helper
    size_t bBytes() const {
        return (size_t)b_size_*b_size_*num_atoms_*num_pj_*sizeof(Real);
    }

    // Read-only inputs
    const Real* model_in_;
    const int*  atoms_;
    const Real* mea_proj_;
    const Real* b_;
    const Real* h_;
    const Real* Rs_;

    // Persistent GPU buffers
    Real* model_;
    Real* proj_;
    Real* proj_diff_;
    Real* Grad_x_;
    Real* Grad_y_;
    Real* Grad_z_;
    Real* Grad_model_;

    // FFT filter
    // std::unique_ptr<FFTFilterR2C> filter_;
    FFTFilterR2C* filter_;

    // Geometry
    dim3 grid_fwd_, block_fwd_;
    dim3 grid_diff_, block1D_;

    // Parameters
    int num_atoms_, num_pj_, num_pts_;
    int b_size_, N1_, N2_, halfWidth_;
    int num_threads_;
    float dt_, scale_;
};