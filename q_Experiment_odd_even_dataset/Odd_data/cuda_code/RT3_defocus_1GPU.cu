#include <iostream>
#include "mex.h"
#include "omp.h"
#include <cmath>
#include "gpu/mxGPUArray.h"

#include <cuda_runtime.h>
#include "utility.cu"

// #include "helper_functions.h"
// #include "helper_cuda.h"
// #include "convolutionSeparable_common.h"


using namespace std;

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

__device__ static float atomicMax(float* address, float val)
{
    int* address_as_i = (int*) address;
    int old = *address_as_i, assumed;
    do {
        assumed = old;
        old = ::atomicCAS(address_as_i, assumed,
            __float_as_int(::fmaxf(val, __int_as_float(assumed))));
    } while (assumed != old);
    return __int_as_float(old);
}

__device__ static float atomicMin(float* address, float val)
{
    int* address_as_i = (int*) address;
    int old = *address_as_i, assumed;
    do {
        assumed = old;
        old = ::atomicCAS(address_as_i, assumed,
            __float_as_int(::fminf(val, __int_as_float(assumed))));
    } while (assumed != old);
    return __int_as_float(old);
}

__device__ __forceinline__ float atomicMinFloat (float * addr, float value) {
        float old;
        old = (value >= 0) ? __int_as_float(atomicMin((int *)addr, __float_as_int(value))) :
             __uint_as_float(atomicMax((unsigned int *)addr, __float_as_uint(value)));

        return old;
}

__device__ __forceinline__ float atomicMaxFloat (float * addr, float value) {
    float old;
    old = (value >= 0) ? __int_as_float(atomicMax((int *)addr, __float_as_int(value))) :
         __uint_as_float(atomicMin((unsigned int *)addr, __float_as_uint(value)));

    return old;
}



template <typename T>
void __global__ updateRec(T*rec_new, int dimx, int dimy, T scale_factor){
    //long long const i = blockDim.x * blockIdx.x + threadIdx.x;
    int const x = blockIdx.x * blockDim.x + threadIdx.x;;
    int const y = blockIdx.y;
    int const z = blockIdx.z;
    if (x<dimx ) {
        long long i = (long long)z*dimx*dimy + y*dimx +  x;
        rec_new[i] = max( 0.0, rec_new[i] ) * scale_factor;
    }
}

template <typename T>
void __global__ compute_xyz_shift( const T*Matrix, const T* shift,  T*x_shift, T*y_shift, T*z_shift, int Num_pjs){
    int const i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<Num_pjs ) {
        int index = 9*i;
        for (int j=0; j<4; j++){
            x_shift[4*i+j] = Matrix[index+0]*shift[2*j] + Matrix[index+3]*0.0 + Matrix[index+6]*shift[2*j+1] ;
            y_shift[4*i+j] = Matrix[index+1]*shift[2*j] + Matrix[index+4]*0.0 + Matrix[index+7]*shift[2*j+1] ;
            z_shift[4*i+j] = Matrix[index+2]*shift[2*j] + Matrix[index+5]*0.0 + Matrix[index+8]*shift[2*j+1] ;
        }
    }   
    
}

template <typename T>
void __global__ R1norm(const T *d_vec, double* R1, int N){
    long long const i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N ) {        
        atomicAdd( R1 , (double)abs(d_vec[i]) );
    }
}

template <typename T>
void __global__ R1norm(const T *d_vec, T* R1, int N){
    long long const i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N ) {        
        atomicAdd( R1 , abs(d_vec[i]) );
    }
}

template <typename T>
void __global__ computeCoords(T* coords, const int dimx, const int dimy, const int ncx, const int ncy, int ncz, long long N, long long starting_point){
    long long i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N ) { 
        i+=starting_point;
        coords[3*i]   =  int(           i%dimx ) - ncx  + 1 ;
        //coords[3*i+1] =  int(  ( i%(dimx*dimy) ) /dimx ) - ncy + 1;
        coords[3*i+1] =  ( int( i/dimx ) ) % dimy - ncy + 1;
        coords[3*i+2] =  int(    i/(dimx*dimy) ) - ncz + 1 ;
    }    
}
//xx_h = mod( ((0:692)'),7)-3;
//yy_h = floor( mod( (0:692)', 7*9) /7)-4; yy_h = mod( floor( (0:692)' /7), 9 )-4;
//zz_h = floor( (0:692)'/(7*9)) - 5;

template <typename T>
void __global__ setValue(T*residual, double val, int N){
    long long const i = blockDim.x * blockIdx.x + threadIdx.x;
    //T o_ratio_inv = 1.0/o_ratio;
    if (i<N) {
        residual[i] = val;
    }
}

//static const int blockSize = 1024;
//static const int gridSize = 24; 



template <typename T>
void __global__ computeResidual(T*residual, const T* projections, const T scale, long long N){
    long long const i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N) {
        residual[i] = residual[i]*scale - projections[i];
    }
}

template <typename T>
void __global__ compute_z_coord(T*z_coord,const T*Matrix, const int nrows, const int ncols, const int dimz,const T*nc, const int o_ratio,const T*x_s,const T*y_s,const T*z_s){
    //long long const i = blockDim.x * blockIdx.x + threadIdx.x;    
    int const x = blockIdx.x * blockDim.x + threadIdx.x;;
    int const y = blockIdx.y;
    int const z = blockIdx.z;
    
    //int origin_offset = 1;
    //long s ;  
    //#pragma omp parallel for default(shared) private(i,s) schedule(static)
    if (x<nrows){
        //const T coord_x = int(           i%nrows ) - nc[0]  + 1 ;
        //const T coord_y = ( int( i/nrows ) ) % ncols - nc[1] + 1;
        const T coord_x = x - nc[0] + 1;;
        const T coord_y = y - nc[1] + 1;
        const T coord_z = z - nc[2] + 1 ;

        //long index = i*3;
        //const T x_i = Matrix[0]*coord_x + Matrix[3]*coord_y + Matrix[6]*coord_z + nc[0];
        //const T y_i = Matrix[1]*coord_x + Matrix[4]*coord_y + Matrix[7]*coord_z + nc[1];
        const T z_c = Matrix[2]*coord_x + Matrix[5]*coord_y + Matrix[8]*coord_z;

        //T x_is = x_i - origin_offset;
        //T y_is = y_i - origin_offset;

        atomicMax( &z_coord[x + y*nrows]              , z_c );
        atomicMin( &z_coord[x + y*nrows + nrows*ncols], z_c );
        
    }
}

template <typename T>
void __global__ sum_slices( T* residual, const T* slices, int dimx, int dimy , int num_slices){
    int const x = blockIdx.x * blockDim.x + threadIdx.x;;
    int const y = blockIdx.y;    
    if(x<dimx){
        double val=0;
        long long N = (long long int) dimx * dimy;
        for(int s=0; s<num_slices; s++){
            val += double( slices[x + y*dimx + s*N] );
        }
        residual[x + y*dimx] = val;
    }    
}

template <typename T>
void __global__ RadonTF(const T*data,const T*Matrix, const int nrows, const int ncols, const int dimz, const T*nc, const int o_ratio, const T*x_s,const T*y_s,
T*slices, int k_min, int k_max, const T defocus_step){
    //long long const i = blockDim.x * blockIdx.x + threadIdx.x;    
    int const x = blockIdx.x * blockDim.x + threadIdx.x;;
    int const y = blockIdx.y;
    int const z = blockIdx.z;
    
    int origin_offset = 1;
    long s ;  
    //#pragma omp parallel for default(shared) private(i,s) schedule(static)
    if (x<nrows){
        const T & data_i = data[z*ncols*nrows + y*nrows + x];
        //const T coord_x = int(           i%nrows ) - nc[0]  + 1 ;
        //const T coord_y = ( int( i/nrows ) ) % ncols - nc[1] + 1;
        const T coord_x = x - nc[0] + 1;;
        const T coord_y = y - nc[1] + 1;
        const T coord_z = z - nc[2] + 1 ;

        //long index = i*3;
        const T x_i = Matrix[0]*coord_x + Matrix[3]*coord_y + Matrix[6]*coord_z + nc[0];
        const T y_i = Matrix[1]*coord_x + Matrix[4]*coord_y + Matrix[7]*coord_z + nc[1];
        const T z_c = Matrix[2]*coord_x + Matrix[5]*coord_y + Matrix[8]*coord_z;
        int k_slice = roundf( z_c/defocus_step );
        k_slice = max( min(k_slice, k_max), k_min);
        T * result = slices + (k_slice-k_min)* nrows*ncols ;

        for (s=0; s<o_ratio; s++){
            //for (i = 0; i < N; ++i) {
            T x_is = x_i + x_s[s] - origin_offset;
            T y_is = y_i + y_s[s] - origin_offset;

            // get coordinates of bounding grid locations
            long long x_1 = ( long long) floor(x_is) ;
            long long x_2 = x_1 + 1;
            long long y_1 = ( long long) floor(y_is) ;
            long long y_2 = y_1 + 1;
            
            if (x_1>=-1 && x_2<=nrows  &&  y_1>=-1 && y_2<=ncols ){ 
                T w_x1 = x_2 - x_is ;
                T w_x2 = 1   - w_x1;
                T w_y1 = y_2 - y_is ;
                T w_y2 = 1   - w_y1;            
                if (x_1==-1){
                    if (y_1==-1){
                        atomicAdd( &result[x_2 + y_2*nrows] , w_x2*w_y2 * data_i);
                    }
                    else if(y_2==ncols){
                        atomicAdd( &result[x_2 + y_1*nrows] , w_x2*w_y1 * data_i);
                    }
                    else{
                        atomicAdd( &result[x_2 + y_1*nrows] , w_x2*w_y1 * data_i);
                        atomicAdd( &result[x_2 + y_2*nrows] , w_x2*w_y2 * data_i);                    
                    }
                }
                else if (x_2==nrows){
                    if (y_1==-1){
                        atomicAdd( &result[x_1 + y_2*nrows] , w_x1*w_y2 * data_i);
                    }
                    else if(y_2==ncols){
                        atomicAdd( &result[x_1 + y_1*nrows] , w_x1*w_y1 * data_i);
                    }
                    else{
                        atomicAdd( &result[x_1 + y_1*nrows] , w_x1*w_y1 * data_i);
                        atomicAdd( &result[x_1 + y_2*nrows] , w_x1*w_y2 * data_i);                  
                    } 
                }
                else{
                    if (y_1==-1){
                        atomicAdd( &result[x_1 + y_2*nrows] , w_x1*w_y2 * data_i);
                        atomicAdd( &result[x_2 + y_2*nrows] , w_x2*w_y2 * data_i);
                    }
                    else if(y_2==ncols){
                        atomicAdd( &result[x_1 + y_1*nrows] , w_x1*w_y1 * data_i);
                        atomicAdd( &result[x_2 + y_1*nrows] , w_x2*w_y1 * data_i);
                    }
                    else{
                        atomicAdd( &result[x_1 + y_1*nrows] , w_x1*w_y1 * data_i);
                        atomicAdd( &result[x_1 + y_2*nrows] , w_x1*w_y2 * data_i);
                        atomicAdd( &result[x_2 + y_1*nrows] , w_x2*w_y1 * data_i);
                        atomicAdd( &result[x_2 + y_2*nrows] , w_x2*w_y2 * data_i);                  
                    }                               
                }
            }
        }
    }
}



template <typename T>
void __global__ RadonTpose_updateRec(const T* Matrix, const int nrows, const int ncols, const T* nc, const T* slices, 
const int o_ratio, const T*x_s,const T*y_s, T* Rec, float dt, long long N, int k_min, int k_max, const T defocus_step){
    //long long const i = blockDim.x * blockIdx.x + threadIdx.x;
    int const x = blockIdx.x * blockDim.x + threadIdx.x;
    int const y = blockIdx.y;
    int const z = blockIdx.z;
    int origin_offset = 1;
    long s;
    //#pragma omp parallel for default(shared) private(s) schedule(static)  
    if( x < nrows ) {
        //const T coord_x = int(           i%nrows )   - nc[0] + 1 ;
        //const T coord_y = ( int( i/nrows ) ) % ncols - nc[1] + 1;
        //const T coord_z =  int(    i/(nrows*ncols) ) - nc[2] + 1 ;
        const T coord_x = x - nc[0] + 1;;
        const T coord_y = y - nc[1] + 1;
        const T coord_z = z - nc[2] + 1 ;
        long long i = (long long)z*ncols*nrows + y*nrows + x;

        //long index = i*3;
        const T x0  = Matrix[0]*coord_x + Matrix[3]*coord_y + Matrix[6]*coord_z + nc[0];
        const T y0  = Matrix[1]*coord_x + Matrix[4]*coord_y + Matrix[7]*coord_z + nc[1];
        const T z_c = Matrix[2]*coord_x + Matrix[5]*coord_y + Matrix[8]*coord_z;
        int k_slice = roundf( z_c/defocus_step );
        k_slice = max( min(k_slice, k_max), k_min);
        const T * data = slices + (k_slice-k_min)* nrows*ncols ;

        for (s=0; s<o_ratio; s++){
            T x_i = x0 + x_s[s];
            T y_i = y0 + y_s[s];
            // get coordinates of bounding grid locations
            long long x_1 = ( long long) floor(x_i) - origin_offset;
            long long x_2 = x_1 + 1;
            long long y_1 = ( long long) floor(y_i) - origin_offset;
            long long y_2 = y_1 + 1;
            
            // handle special case where x/y is the last element
            if ( (x_i - origin_offset) == (nrows-1) )   { x_2 -= 1; x_1 -= 1;}
            if ( (y_i - origin_offset) == (ncols-1) )   { y_2 -= 1; y_1 -= 1;}
            
            // return 0 for target values that are out of bounds
            if (x_1 < 0 | x_2 > (nrows - 1) |  y_1 < 0 | y_2 > (ncols - 1)){
                //result[i] = 0;
            }
            else {
                // get the array values
                const T& f_11 = data[x_1 + y_1*nrows];
                const T& f_12 = data[x_1 + y_2*nrows];
                const T& f_21 = data[x_2 + y_1*nrows];
                const T& f_22 = data[x_2 + y_2*nrows];
                
                // compute weights
                T w_x1 = x_2 - (x_i - origin_offset);
                T w_x2 = (x_i - origin_offset) - x_1;
                T w_y1 = y_2 - (y_i - origin_offset);
                T w_y2 = (y_i - origin_offset) - y_1;
                
                T a,b;
                a = f_11 * w_x1 + f_21 * w_x2;
                b = f_12 * w_x1 + f_22 * w_x2;
                Rec[i] -= dt*(a * w_y1 + b * w_y2);
            }
        }
    }
}

template <typename T>
void __global__ compute_kR_rad2(T*kR_rad2, T dkx, const int nr, const int c_nr, const long long nr2){
    long long const ii = blockIdx.x * blockDim.x + threadIdx.x;
    if(ii<nr2){
        const int x = ii/nr - c_nr;
        const int y = ii%nr - c_nr;
        kR_rad2[ii] = (x*x + y*y)*dkx*dkx;
    }
}

template <typename T>
void __global__ compute_ape(bool*ape, const T*kR_rad2, T alpha2, const long long nr2){
    long long const ii = blockIdx.x * blockDim.x + threadIdx.x;
    if(ii<nr2){
        ape[ii] = kR_rad2[ii]<=alpha2;
    }
}

template <typename T>
void __global__ avg_y(const T* d_slices, T*d_slices_tmp, const int dimx, const int dimy, const int num_slices){
    int i =  blockIdx.x * blockDim.x + threadIdx.x;
    int j =  blockIdx.y * blockDim.y + threadIdx.y;
    if(i<dimx && j<dimy){
        d_slices_tmp += (i+dimx*j);
        d_slices     += (i+dimx*j);
        if(j==0){
            for(int k=0; k<num_slices; k++){
                d_slices_tmp[0] = 0.8571429f* d_slices[0] + 0.1428571f* d_slices[dimx];
                d_slices_tmp += dimx*dimy;
                d_slices     += dimx*dimy;
            }      
        }
        else if(j==dimy-1){
            for(int k=0; k<num_slices; k++){
                d_slices_tmp[0] = 0.8571429f* d_slices[0] + 0.1428571f* d_slices[-dimx];
                d_slices_tmp += dimx*dimy;
                d_slices     += dimx*dimy;
            }        
        }
        else{
            for(int k=0; k<num_slices; k++){
                d_slices_tmp[0] =  0.75f* d_slices[0] + 0.125f* d_slices[dimx] + 0.125* d_slices[-dimx];
                d_slices_tmp += dimx*dimy;
                d_slices     += dimx*dimy;
            }        
        }
    }
}


void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]) {
    int const threadsPerBlock = 256;
    int blocksPerGridPrj; // blocksPerGridRec;

    /* Initialize the MathWorks GPU API. */
    mxInitGPU();    
    cudaError_t err = cudaSuccess;
    err = cudaSetDevice(0);            // Set device 0 as current    
    if(err!=cudaSuccess){
        printf("cuda fail to set\n");
    }

    int deviceCount;
    cudaGetDeviceCount(&deviceCount);
    int device;
    for (device = 0; device < deviceCount; ++device) {
        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, device);
        printf("Device %d has compute capability %d.%d.\n",
        device, deviceProp.major, deviceProp.minor);
    }

    char const * const errId = "parallel:gpu:mexGPUExample:InvalidInput";
    char const * const errMsg = "Invalid input to MEX file.";    
    if ((nrhs!=7)&&(nrhs!=8) ){ //|| (nlhs!=1)  ) {
        cout << "error: number of parameter not matched"<<endl;
        mexErrMsgIdAndTxt(errId, errMsg); //!(mxIsGPUArray(prhs[0]))
    }
    /*
    * 0: projections
    * 1: rotmat
    * 2: dimz
    * 3: iterations
    * 4: step_size
    * 5: positivity?
    * 6: defocus_step
    * 7: Defocus_param
    * 8: semi_angle
    */

    /*
    float *h_Kernel, *h_Input, *h_Buffer, *h_OutputCPU, *h_OutputGPU;
    float *d_Input, *d_Output, *d_Buffer, * d_kernel;

    const int imageW = 3072;
    const int imageH = 3072;
    const int iters = 16;

    StopWatchInterface *hTimer = NULL;
    sdkCreateTimer(&hTimer);

    printf("Image Width x Height = %i x %i\n\n", imageW, imageH);
    printf("Allocating and initializing host arrays...\n");
    h_Kernel = (float *)malloc(KERNEL_LENGTH * sizeof(float));
    h_Input = (float *)malloc(imageW * imageH * sizeof(float));
    h_Buffer = (float *)malloc(imageW * imageH * sizeof(float));
    h_OutputCPU = (float *)malloc(imageW * imageH * sizeof(float));
    h_OutputGPU = (float *)malloc(imageW * imageH * sizeof(float));
    srand(200);

    for (unsigned int i = 0; i < KERNEL_LENGTH; i++) {
        h_Kernel[i] = (float)(rand() % 16);
    }

    for (unsigned i = 0; i < imageW * imageH; i++) {
        h_Input[i] = (float)(rand() % 16);
    }

    printf("Allocating and initializing CUDA arrays...\n");
    checkCudaErrors(
        cudaMalloc((void **)&d_Input, imageW * imageH * sizeof(float)));
    checkCudaErrors(
        cudaMalloc((void **)&d_Output, imageW * imageH * sizeof(float)));
    checkCudaErrors(
        cudaMalloc((void **)&d_Buffer, imageW * imageH * sizeof(float)));
    checkCudaErrors(
        cudaMalloc((void **)&d_kernel, KERNEL_LENGTH* sizeof(float)));
    checkCudaErrors(
        cudaMemcpy( d_kernel, h_Kernel, KERNEL_LENGTH* sizeof(float), cudaMemcpyHostToDevice) );

    setConvolutionKernel(h_Kernel);
    checkCudaErrors(cudaMemcpy(d_Input, h_Input, imageW * imageH * sizeof(float),
                                cudaMemcpyHostToDevice));

    printf("Running GPU convolution (%u identical iterations)...\n\n",
            iters);

    //mexMatlabcall('probefilter')
    for (int i = -1; i < iters; i++) {
        // i == -1 -- warmup iteration
        if (i == 0) {
        checkCudaErrors(cudaDeviceSynchronize());
        sdkResetTimer(&hTimer);
        sdkStartTimer(&hTimer);
        }
        convolutionRowsGPU(d_Buffer, d_Input, imageW, imageH, d_kernel);
        convolutionColumnsGPU(d_Output, d_Buffer, imageW, imageH);
    }

    checkCudaErrors(cudaDeviceSynchronize());
    sdkStopTimer(&hTimer);
    double gpuTime = 0.001 * sdkGetTimerValue(&hTimer) / (double)iters;
    printf(
        "convolutionSeparable, Throughput = %.4f MPixels/sec, Time = %.5f s, "
        "Size = %u Pixels, NumDevsUsed = %i, Workgroup = %u\n",
        (1.0e-6 * (double)(imageW * imageH) / gpuTime), gpuTime,
        (imageW * imageH), 1, 0);

    printf("\nReading back GPU results...\n\n");
    checkCudaErrors(cudaMemcpy(h_OutputGPU, d_Output,
                                imageW * imageH * sizeof(float),
                                cudaMemcpyDeviceToHost));
    */

    // h = createProbeKernel(V, pixelsize, nr, alpha, df_series);
    
    
    mxArray *lhs[1], *rhs[5];
    // double Voltage = 300e3;
    // double pixelsize = 0.3079;
    // double nr=200;
    // double semi_ang = 22e-3;
    // double defocusstep = 20;

    //mxArray *rhs_4 = mxCreateDoubleScalar(df_series);
    //const mwSize res_Size[] = {1};
    //const mxGPUArray * res = mxGPUCreateGPUArray(mwSize(1), res_Size, mxDOUBLE_CLASS, mxREAL, MX_GPU_INITIALIZE_VALUES );
    //rhs[4] = mxGPUCreateMxArrayOnGPU(res);
    //double *res_pr = mxGetDoubles( (mxArray*)res);
    //lhs[0] = mxCreateDoubleMatrix(m, n , mxREAL);
    //lhs[0] = mxCreateNumericArray(3, recSize, mxDOUBLE_CLASS, mxREAL);
    //double *lhs_0 = mxGetPr(lhs[0]);



    

    // read inputs
    const float  * projections  = mxGetSingles(prhs[0]);
    const float  * Matrix       = mxGetSingles(prhs[1]);
    const double * dimzPtr      = mxGetPr(prhs[2]);  const mwSize dimz = mwSize( dimzPtr[0] ); 

    const double * GD_info = mxGetPr(prhs[3]);  
    int iterations         = int( GD_info[0] ); 
    double step_size       = double( GD_info[1] ); 

    const double * constraints = mxGetPr(prhs[4]);  //mxGetLogicals(prhs[4]); 
    const bool positivity      = bool( constraints[0] ); 
    const bool is_y_avg        = constraints[1] != 0.0f;
    const float l2_reg         = constraints[2];
    float scale_factor = 1.0f / ( 1.0f + step_size*l2_reg );    
    
    //defocus_info = [Voltage, pixelsize, nr, semi_angle, defocus_step];
    const double * defocus_info  = mxGetPr(prhs[5]);  
    double Voltage       = defocus_info[0]; 
    double pixelsize     = defocus_info[1]; 
    unsigned long nr     = (unsigned long) ( defocus_info[2] ); 
    double semi_angle    = defocus_info[3]; 
    float defocus_step   = float( defocus_info[4] ); 
    float defocus_scale  = float( defocus_info[5] );

    const float  * defocus_param = mxGetSingles(prhs[6]);


    double alpha = semi_angle*1e3;
    double df_series = abs(2*defocus_step-0)*pixelsize;
    rhs[0] = mxCreateDoubleScalar(Voltage);
    rhs[1] = mxCreateDoubleScalar(pixelsize);
    rhs[2] = mxCreateDoubleScalar(nr);
    rhs[3] = mxCreateDoubleScalar(alpha);
    rhs[4] = mxCreateDoubleScalar(df_series);

    // mexCallMATLAB(1, lhs, 5, rhs, "createProbeKernel");
    // plhs[2] = lhs[0];

    const size_t o_ratio=4;
    const mwSize* projsSize   = (mxGetDimensions(prhs[0]));
    const mwSize dimx    = projsSize[0];
    const mwSize dimy    = projsSize[1];
    const mwSize Num_pjs = projsSize[2];

    const long long nrow_cols  = dimx*dimy;
    const long long nPjsPoints = dimx*dimy*Num_pjs;
    const long long recPoints = dimx*dimy*dimz;
    const mwSize recSize[] = {dimx,dimy,dimz};
    //cout << dimx <<", " << dimy <<", " << dimz <<", "  << endl;

    const mwSize* dims_Mat = (mxGetDimensions(prhs[1]));
    const mwSize Num_pjs2   = dims_Mat[2];
    const mwSize R1 = dims_Mat[0];
    const mwSize R2 = dims_Mat[1];

    const mwSize* dims_defocus_param = (mxGetDimensions(prhs[6]));
    const mwSize Num_pjs3   = dims_defocus_param[0];

    /*
    const mwSize* o_ratio_dim = mxGPUGetDimensions(xs2);
    const size_t o_ratio      = o_ratio_dim[0];
    const size_t Num_pjs2      = o_ratio_dim[1];                
    */
    if(Num_pjs2!=Num_pjs || Num_pjs3!=Num_pjs || R1!=3 || R2!=3 ){
        cout << "error: dimension not matched"<<endl;
        mexErrMsgIdAndTxt(errId, errMsg);
    }
    if(mxGetClassID(prhs[0]) == mxDOUBLE_CLASS || mxGetClassID(prhs[1]) == mxDOUBLE_CLASS ){     
        printf("can't work with double\n");
        return;
    }    

    const mwSize ncx = int (floor(dimx/2.0)+1);
    const mwSize ncy = int (floor(dimy/2.0)+1);
    const mwSize ncz = int (floor(dimz/2.0)+1);

    const double dt =  (step_size/Num_pjs/dimz/o_ratio);
    //mexPrintf("%d\n",npoints);   
    
    // copy projections to GPU
    float * d_projections;
    cudaMalloc( &d_projections, nPjsPoints*sizeof(float) );
    cudaMemcpy(  d_projections, projections, nPjsPoints*sizeof(float), cudaMemcpyHostToDevice );

    // copy rotation matrix to GPU
    float * d_Matrix;
    cudaMalloc( &d_Matrix, 9*Num_pjs*sizeof(float) );
    cudaMemcpy( d_Matrix, Matrix, 9*Num_pjs*sizeof(float), cudaMemcpyHostToDevice );


    // plhs[0] is the returning reconstruction
    plhs[0] = mxCreateNumericArray(3, recSize, mxSINGLE_CLASS, mxREAL);
    float * rec   = (float*)mxGetSingles(plhs[0]);

    // plhs[1] is the returning calculated projection
    plhs[1] = mxCreateNumericArray(3, projsSize, mxSINGLE_CLASS, mxREAL);
    float* cal_proj   = (float*)mxGetSingles(plhs[1]);

    // create reconstruciton on GPU
    float * d_rec;
    cudaMalloc( &d_rec,      recPoints*sizeof(float) );
    cudaMemset( d_rec, 0.0f, recPoints*sizeof(float));

    if(nrhs==8){
        const mwSize* recSize2   = mxGetDimensions(prhs[7]);
        const mwSize dimx2   = recSize2[0];
        const mwSize dimy2   = recSize2[1];
        const mwSize dimz2   = recSize2[2];
        if( dimx2!=dimx || dimy2!=dimy || dimz2!=dimz ){
            cout << "error: reconstruction shape not matched"<<endl;
            mexErrMsgIdAndTxt(errId, errMsg);
        }
        float * rec_ori = (float*)mxGetSingles(prhs[7]);        
        //memcpy( h_rec, rec_ori, recPoints*sizeof(float) );
        cudaMemcpy( d_rec, rec_ori, recPoints*sizeof(float), cudaMemcpyHostToDevice);
    }     
    else{
        cudaMemset(d_rec, 0, recPoints*sizeof(float));
    }
    
    // create residual on GPU
    float * d_residual;
    cudaMalloc( &d_residual, nPjsPoints*sizeof(float) );
    /*
    mxGPUArray * residual = mxGPUCreateGPUArray(mwSize(3),
                        pjsSize, mxGPUGetClassID(projections),
                        mxREAL, MX_GPU_INITIALIZE_VALUES );
    */
    
    
    //blocksPerGrid    = (recPoints + threadsPerBlock - 1) / threadsPerBlock;
    blocksPerGridPrj = (nPjsPoints + threadsPerBlock - 1) / threadsPerBlock;
    //blocksPerGridRec = (recPoints + threadsPerBlock - 1) / threadsPerBlock;
    dim3 blocksPerGridRec2( (dimx + threadsPerBlock - 1) / threadsPerBlock , dimy, dimz );

    //mxGPUArray const *support     = mxGPUCreateFromMxArray(prhs[7]);
    //bool const *d_support = (bool const *) (mxGPUGetDataReadOnly(support));            


    //float const *d_projections = (float const *) (mxGPUGetDataReadOnly(projections));
    //float const *d_Matrix      = (float const *) (mxGPUGetDataReadOnly(Matrix));
    //float const *d_Coord       = (float const *) (mxGPUGetDataReadOnly(Coord));
    //float const *d_xs         = (float const *) (mxGPUGetDataReadOnly(xs2));
    //float const *d_ys         = (float const *) (mxGPUGetDataReadOnly(ys2));
    //float *d_residual = (float *) (mxGPUGetData(residual));


    // compute rotated shift
    float shift[]  = {0.25,0.25, 0.25,-0.25, -0.25,0.25, -0.25,-0.25};
    float *shift_ptr;    
    cudaMalloc( (void**) &shift_ptr, 8*sizeof(float) );
    cudaMemcpy( shift_ptr, shift, 8*sizeof(float), cudaMemcpyHostToDevice);
    float * x_shift, *y_shift, * z_shift;
    cudaMalloc( (void**) &x_shift, 4*Num_pjs*sizeof(float) );
    cudaMalloc( (void**) &y_shift, 4*Num_pjs*sizeof(float) );
    cudaMalloc( (void**) &z_shift, 4*Num_pjs*sizeof(float) );
    compute_xyz_shift<<<2, threadsPerBlock>>>( d_Matrix, shift_ptr, x_shift, y_shift, z_shift, Num_pjs );
    float const *d_xs2         = (float  *) x_shift;
    float const *d_ys2         = (float  *) y_shift;
    float const *d_zs2         = (float  *) z_shift;

    // compute cartesian coordinates
    //computeCoords<<<blocksPerGridRec,threadsPerBlock>>>(d_Coord, dimx, dimy, ncx,ncy, ncz, recPoints, 0);

    // compute nc = [ncx,ncy,ncz]
    //const float nc_cpu[]  = { ncx,ncy,ncz}; 
    const float nc_cpu[]  = { float(floor(dimx/2.0)+1), float(floor(dimy/2.0)+1), float(floor(dimz/2.0)+1)};     
    float * d_nc;
    cudaMalloc( (void**)&d_nc, 3*sizeof(float) );
    cudaMemcpy( d_nc, nc_cpu,  3*sizeof(float), cudaMemcpyHostToDevice ); 

    // compute norm of projection
    /*double pj_norm,  *pj_norm_pt;
    cudaMalloc( (void**)&pj_norm_pt, sizeof(double) );
    cudaMemset(pj_norm_pt, 0, 1*sizeof(double));   //setValue<<<1,1>>>(pj_pt, 0.0, 1);
    R1norm<<<blocksPerGridPrj, threadsPerBlock>>>(d_projections,  pj_norm_pt, nPjsPoints);
    cudaMemcpy( &pj_norm ,  pj_norm_pt, sizeof(double), cudaMemcpyDeviceToHost ) ;
    cudaFree( pj_norm_pt);
    */
    
    double pj_norm = sumArray(d_projections,nPjsPoints);

    // create temporary z_coord
    float * d_z_coord;
    //float * z_coord_1 = new float;
    //float * z_coord_2 = new float;
    cudaMalloc( &d_z_coord,       2*((long long) dimx)*dimy*sizeof(float) );
    cudaMemset(  d_z_coord, 0.0f, 2*((long long) dimx)*dimy*sizeof(float) );
    //cudaMemcpy(  d_projections, projections, nPjsPoints*sizeof(float), cudaMemcpyHostToDevice );
    int * pj_slices = new int[Num_pjs*2]{};
    for(int i=0; i<Num_pjs; i++){
        //cudaDeviceSynchronize();
        compute_z_coord<<<blocksPerGridRec2, threadsPerBlock>>>(d_z_coord, d_Matrix + i*9, dimx,dimy,dimz,d_nc,
            o_ratio, d_xs2+i*o_ratio,d_ys2+i*o_ratio,d_zs2+i*o_ratio);        
        double z_max = maxArray(d_z_coord , ( (long long) dimx )*dimy);
        double z_min = minArray(d_z_coord + dimx*dimy , ( (long long) dimx )*dimy);   
        //cudaMemcpy( z_coord_1, d_z_coord, sizeof(float), cudaMemcpyDeviceToHost);    
        //cudaMemcpy( z_coord_2, d_z_coord + dimx*dimy, sizeof(float), cudaMemcpyDeviceToHost);    
        int k_min = round(z_min/defocus_step);
        int k_max = round(z_max/defocus_step);
        //int num_k = k_max-k_min + 1;
        pj_slices[2*i]   = k_min;
        pj_slices[2*i+1] = k_max;
        cudaMemset(d_z_coord, 0.0f, 2*dimx*dimy*sizeof(float) );
        //cout << i <<". z_min=" << z_min << ", z_max=" << z_max << ", num_k="<<num_k <<endl;
    }
    cudaFree(d_z_coord);    
    dim3 blocksPerGridSlice( (dimx + threadsPerBlock - 1) / threadsPerBlock , dimy );

    //cout << "dimx=" <<dimx <<", dimy=" << dimy << endl;

    //defocus setting
    const double hh = 6.62e-34;
    const double ee = 1.6e-19;
    const double mm = 9.1e-31;
    const double cc = 3e8;
    const double lambda = 1e10*hh/sqrt(2*mm*ee*Voltage) / sqrt(1 + ee*Voltage/ (2.0*mm*cc*cc) );
    double kmax = 1.0e3/pixelsize*lambda/2;
    double dkx = 2*kmax/nr;

    // call Matlab function myFrequencyProbeFilter
    // call FrequencyProbeFilter_forward(A, kR_rad2, ape, df_series, lambda , hsize)
    mxArray *LHS[1], *RHS[7];

    const mwSize nrSize[] = {nr+1,nr+1};
    RHS[1] = mxCreateNumericArray(2, nrSize, mxSINGLE_CLASS, mxREAL);
    float * kR_rad2   = (float*)mxGetSingles(RHS[1]);
    float * d_kR_rad2;
    long long nr2 = (long long)(nr+1)*(nr+1);
    const int c_nr = (nr+1)/2;
    cudaMalloc( &d_kR_rad2, nr2*sizeof(float) );
    cudaMemset(  d_kR_rad2, 0, nr2*sizeof(float) );
    compute_kR_rad2<<< (nr2+threadsPerBlock-1)/threadsPerBlock, threadsPerBlock >>> (d_kR_rad2, (float)dkx, nr+1, c_nr, nr2);
    cudaMemcpy( kR_rad2, d_kR_rad2, nr2*sizeof(float), cudaMemcpyDeviceToHost);

    RHS[2] = mxCreateLogicalMatrix(nr+1, nr+1);
    bool * ape   = (bool*)mxGetLogicals(RHS[2]);
    float alpha2 = alpha*alpha;
    bool * d_ape;
    cudaMalloc( &d_ape, nr2*sizeof(bool) );
    cudaMemset(  d_ape, 0, nr2*sizeof(bool) );
    compute_ape <<<(nr2+threadsPerBlock-1)/threadsPerBlock, threadsPerBlock >>> (d_ape, d_kR_rad2, alpha2, nr2);
    cudaMemcpy( ape, d_ape, nr2*sizeof(bool), cudaMemcpyDeviceToHost);
    
    RHS[4] = mxCreateDoubleScalar(lambda);

    mwSize hsize_size[] = {1,2};
    RHS[5] = mxCreateNumericArray(2, hsize_size, mxSINGLE_CLASS, mxREAL);
    float * hsize   = (float*)mxGetSingles(RHS[5]);
    hsize[0] = nr;
    hsize[1] = nr;
    
    //myFrequencyProbeFilter(A, kR_rad2, ape, df_values, lambda, hsize)
    //FrequencyProbeFilter_forward(A, kR_rad2, ape, df_series, lambda , hsize)


    // set zero to residual
    //setValue<<<blocksPerGridPrj, threadsPerBlock>>>(d_residual, 0.0, nPjsPoints);        

    // iteration
    for (long iter=0; iter<iterations; iter++){
        // set zero residual
        //setValue<<<blocksPerGridPrj, threadsPerBlock>>>(d_residual, 0.0, nPjsPoints);
        cudaMemset( d_residual, 0, nPjsPoints*sizeof(float) );

        // compute forward projection        
        for (int i = 0; i < Num_pjs; i++){            
            int k_min = pj_slices[2*i];
            int k_max = pj_slices[2*i+1];
            int num_slices = k_max - k_min +1 ;
            
            float * d_slices, *d_slices_tmp;
            cudaMalloc( &d_slices,    num_slices*((long long) dimx)*dimy*sizeof(float) );            
            cudaMemset(  d_slices, 0, num_slices*((long long) dimx)*dimy*sizeof(float) );
            
            //const mwSize size_slices[3]={dimx,dimy,num_slices};
            //mxGPUArray * mex_slices = mxGPUCreateGPUArray(mwSize(3), size_slices, mxSINGLE_CLASS, mxREAL, MX_GPU_INITIALIZE_VALUES );            
            //float *d_slices = (float*) ( mxGPUGetData( mex_slices) );
            //RHS[0] = (mxArray*) mex_slices;
            //cout << "RHS[0]" <<endl;

            // radon transform for each d_slices
            RadonTF<<<blocksPerGridRec2, threadsPerBlock>>>(d_rec, d_Matrix + i*9, dimx,dimy,dimz,d_nc,
            o_ratio,d_xs2+i*o_ratio,d_ys2+i*o_ratio, d_slices, k_min, k_max, defocus_step );
            //d_residual + i*nrow_cols

            if(is_y_avg){
                cudaMalloc( &d_slices_tmp,    num_slices*((long long) dimx)*dimy*sizeof(float) );
                //cudaMemset(  d_slices_tmp, 0, num_slices*((long long) dimx)*dimy*sizeof(float) );
                avg_y<<< dim3( (dimx+15)/16, (dimy+15)/16 ), dim3(16,16) >>>( d_slices, d_slices_tmp, dimx, dimy, num_slices);
                cudaFree(d_slices);
                d_slices = d_slices_tmp;
            }

            //copy d_slices to RHS[0]:h_slices
            const mwSize size_slices[3]={dimx,dimy,num_slices};
            RHS[0] = mxCreateNumericArray(3, size_slices, mxSINGLE_CLASS, mxREAL);
            float* h_slices = (float*)mxGetSingles(RHS[0]);
            
            cudaMemcpy(h_slices, d_slices, dimx*dimy*num_slices*sizeof(float), cudaMemcpyDeviceToHost);
            if (i==-1){
                for(int x=0; x<dimx; x++){
                    for(int y=0; y<dimy; y++){
                        cout<< h_slices[dimx*x+y + (num_slices/2)*dimx*dimy] <<",";
                    }
                    cout <<endl;
                }
            }

            // call Matlab function
            const mwSize size_df_value[1] = {num_slices};
            RHS[3] = mxCreateNumericArray(1, size_df_value, mxSINGLE_CLASS, mxREAL);
            float * df_values   = (float*)mxGetSingles(RHS[3]);
            //cout<< k_min << "," << k_max << endl;
            for(int k=k_min; k<=k_max; k++){
                df_values[k-k_min] = abs(k*defocus_step-defocus_param[i])*pixelsize*defocus_scale;
                //cout<< k <<"." << df_values[k-k_min]<<", ";
            }
            
            //FrequencyProbeFilter_forward(A, kR_rad2, ape, df_series, lambda , hsize)
            mexCallMATLAB(1, LHS, 6, RHS, "FrequencyProbeFilter_forward"); 
            //plhs[1] = LHS[0];

            // get data from Matlab and copy to gpu
            const float  * h_slices_blurry  = mxGetSingles(LHS[0]);
            cudaMemcpy(d_slices, h_slices_blurry, num_slices*dimx*dimy*sizeof(float),  cudaMemcpyHostToDevice);


            /* to do: convolve d_slices with probe_filter
                for(int k=0, k<num_slices; k++){
                    d_slices[:,:,k] = conv2(d_slices[:,:,k], probe_filter) 
                }
            */


            sum_slices<<<blocksPerGridSlice, threadsPerBlock>>>(d_residual + i*nrow_cols, d_slices, dimx,dimy , num_slices);
            cudaFree(d_slices);
            //mxFree(h_slices_blurry);
            mxDestroyArray(RHS[3]);
            mxDestroyArray(RHS[0]);
            mxDestroyArray(LHS[0]);
        }

        // copy calculated projections
        if(iter==iterations-1){
            cudaMemcpy( cal_proj, d_residual, Num_pjs*dimx*dimy*sizeof(float),  cudaMemcpyDeviceToHost);
        }

        // compute residual: = forward projection - measure projections
        computeResidual<<<blocksPerGridPrj, threadsPerBlock>>>(d_residual,d_projections,1.0f/o_ratio, nPjsPoints);

        // compute R1 factor
        if (iter%1 == 0){
            /*
            double res_norm, *res_norm_pt;
            cudaMalloc( (void**)&res_norm_pt, sizeof(double) );                
            cudaMemset(res_norm_pt, 0, 1*sizeof(double)); //setValue<<<1,1>>>(res_norm_pt, 0.0, 1);
            R1norm<<<blocksPerGridPrj, threadsPerBlock>>>(d_residual,    res_norm_pt, nPjsPoints);
            cudaMemcpy( &res_norm, res_norm_pt, sizeof(double), cudaMemcpyDeviceToHost ) ;
            cudaFree( res_norm_pt );  
            */
            float res_norm, *sum_pt;
            cudaMalloc((void**)&sum_pt, sizeof(float)*gridSize);
            sumCommMultiBlock<<<gridSize, blockSize>>>(d_residual, nPjsPoints, sum_pt);
            sumCommMultiBlock<<<1, blockSize>>>(sum_pt, gridSize, sum_pt);            
            //cudaDeviceSynchronize();
            cudaMemcpy(&res_norm, sum_pt, sizeof(float), cudaMemcpyDeviceToHost);
            cudaFree(sum_pt);

            cout << iter+1 << ". R1 = " << res_norm/pj_norm << endl;
        }            

        // back projection
        for (long i = 0; i < Num_pjs; i ++){            
            int k_min = pj_slices[2*i];
            int k_max = pj_slices[2*i+1];
            int num_slices = k_max - k_min +1 ;
            
            float * d_slices;
            
            cudaMalloc( &d_slices,    num_slices*((long long) dimx)*dimy*sizeof(float) );
            //cudaMemset(  d_slices, 0, num_slices*((long long) dimx)*dimy*sizeof(float) );
            // to do: convolution d_slices with probe_filter instead of copy
            for(int k=0; k<num_slices; k++){
                cudaMemcpy(  d_slices + k*dimx*dimy, d_residual + i*nrow_cols, (long long) dimx * dimy * sizeof(float), cudaMemcpyDeviceToDevice );            
            }
            

            const mwSize size_slices[2]={dimx,dimy};
            RHS[0] = mxCreateNumericArray(2, size_slices, mxSINGLE_CLASS, mxREAL);
            float* h_slices = (float*)mxGetSingles(RHS[0]);
            cudaMemcpy(h_slices, d_residual + i*nrow_cols, dimx*dimy*sizeof(float), cudaMemcpyDeviceToHost);

            // if (iter==4 && i==4){
            //     for(int y=0; y<dimy; y++){
            //         cout<< h_slices[dimx*dimy/2 + y] <<",";
            //     }
            //     cout <<endl;
            // }

            const mwSize size_df_value[1] = {num_slices};
            RHS[3] = mxCreateNumericArray(1, size_df_value, mxSINGLE_CLASS, mxREAL);
            float * df_values   = (float*)mxGetSingles(RHS[3]);
            //cout<< k_min << "," << k_max << endl;
            for(int k=k_min; k<=k_max; k++){
                df_values[k-k_min] = abs(k*defocus_step-defocus_param[i])*pixelsize;
                //cout<< k <<"." << df_values[k-k_min]<<", ";
            }

            //FrequencyProbeFilter_backward(A, kR_rad2, ape, df_series, lambda , hsize)
            mexCallMATLAB(1, LHS, 6, RHS, "FrequencyProbeFilter_backward"); 
            
            // get data from Matlab and copy to gpu
            const float  * h_slices_blurry  = mxGetSingles(LHS[0]);
            // if (iter==4 && i==4){
            //     for(int y=0; y<dimy; y++){
            //         cout<< h_slices_blurry[dimx*dimy/2 + y] <<",";
            //     }
            //     cout <<endl;
            // }       

            cudaMemcpy(d_slices, h_slices_blurry, ((long long) num_slices)*dimx*dimy*sizeof(float), cudaMemcpyHostToDevice);

            RadonTpose_updateRec<<<blocksPerGridRec2, threadsPerBlock>>>(d_Matrix + i*9, dimx,dimy, d_nc, d_slices,
            o_ratio, d_xs2+i*o_ratio, d_ys2+i*o_ratio, d_rec, dt, recPoints, k_min, k_max, defocus_step);
            cudaFree(d_slices);
            //mxFree(h_slices_blurry);
            mxDestroyArray(RHS[3]);
            mxDestroyArray(RHS[0]);
            mxDestroyArray(LHS[0]);
        }

        // update rec
        if(positivity){
            updateRec<<<blocksPerGridRec2, threadsPerBlock>>>(d_rec, dimx, dimy, scale_factor);
        }
    }

    /* return result  */
    //plhs[0] = mxGPUCreateMxArrayOnCPU(rec);
    cudaMemcpy( rec, d_rec, recPoints*sizeof(float), cudaMemcpyDeviceToHost);
    //plhs[1] = mxGPUCreateMxArrayOnGPU(residual);

    //mxGPUDestroyGPUArray(projections); 
    //mxGPUDestroyGPUArray(rec);        
    cudaFree(d_rec);
    cudaFree(d_residual);
    cudaFree(d_Matrix);
    cudaFree(d_projections);
    
    //mxGPUDestroyGPUArray(Matrix);
    //mxGPUDestroyGPUArray(nc);
    cudaFree(d_nc);
    //mxGPUDestroyGPUArray(xs2);
    //mxGPUDestroyGPUArray(ys2);
    //mxGPUDestroyGPUArray(residual);
    //mxGPUDestroyGPUArray(support); 
    cudaFree( shift_ptr );
    cudaFree( x_shift );
    cudaFree( y_shift );
    
    //cudaDeviceReset();      
    //mxGPUArray       *coord2 = (mxGPUArray*)d_Coord2;
    // plhs[1] = mxGPUCreateMxArrayOnCPU(coord2);
    
}


/*
// compute diffnorm of xyshift
float x_norm, x_norm2,  *x_pt, *x_pt2;
cudaMalloc( (void**)&x_pt, sizeof(float) );
cudaMalloc( (void**)&x_pt2, sizeof(float) );
cudaMemset(x_pt, 0, sizeof(float));
cudaMemset(x_pt2, 0, sizeof(float));
R1norm<<<2, threadsPerBlock>>>(d_xs,  x_pt, Num_pjs*4);
R1norm<<<2, threadsPerBlock>>>(d_xs2,  x_pt2, Num_pjs*4);
cudaMemcpy( &x_norm ,  x_pt, sizeof(float), cudaMemcpyDeviceToHost ) ;
cudaMemcpy( &x_norm2 ,  x_pt2, sizeof(float), cudaMemcpyDeviceToHost ) ;
cudaFree( x_pt);
cout << "x_norm = " << x_norm << ", x_norm2 = " << x_norm2 <<endl;
*/



















