% compilation:
% mexcuda -lcufft -R2018a cal_grad_all.cu

mexcuda '-L/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/math_libs/12.4/targets/x86_64-linux/lib' -lcufft -R2018a cal_grad_all.cu
