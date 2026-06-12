clear;clc;
%% add the path for the code
addpath([pwd,'/input_data/'])
addpath([pwd,'/output_data/'])
addpath([pwd,'/src/'])
addpath([pwd,'/src_GPU/'])

%% load data
projections = importdata('Projs_Refinement.mat');
angles      = importdata('angles_Refinement.mat');
model       = importdata('atom_Refinement.mat');
atoms       = importdata('label_Refinement.mat');

[N1,N2,Num_pj] = size(projections);
% the cropped bondary size for each atoms
halfWidth = 10;
% the atomic number for different type:
% use 8 for type 1(O), 22 for type 2(Ti), 56 for type 3(Ba)
Z_arr   = [6, 14];  
% indicate the pixel size for measured projections
Res     = 0.19/2;                      

xdata = [];
xdata.Res       = Res;
xdata.Z_arr     = Z_arr;
xdata.halfWidth = halfWidth;
xdata.atoms     = atoms;
xdata.model     = model;
xdata.angles    = angles;

para0 = [1 1.5;  
         4.5 4.5];

%% refinement
model_refined = model;

for jjjj=1:40
    fprintf('iteration num: %d; \n',jjjj);
    
    xdata.model       = model_refined;
    xdata.model_ori   = model_refined;
    xdata.projections = projections;
    
    tic
    [cal_proj_opt,para_opt] = Cal_HBproj_bayes(para0, xdata, projections);
    toc
    
    para0 = para_opt;
    disp('Initial para0:')
    disp(para0)
    
    xdata.projections = [];    
    xdata.step_sz    = 1;
    xdata.iterations = 2;
    
    tic
    [y_pred,para0,errR] = gradient_B_2type_difB(para0, xdata, projections);
    toc
    
    xdata.step_sz    = 2;
    xdata.iterations = 100;
    
    disp('Updated para0 after gradient_B_2type_difB:')
    disp(para0)
    
    tic
    [para_gpu,errR_gpu, model_out_2_cpu, proj_gpu_2_cpu,grad_x_set_2_cpu, grad_y_set_2_cpu, grad_z_set_2_cpu] = ...
    gradient_XYZ_gpu(para0, xdata, projections);
    toc

    model_refined = model_out_2_cpu;
    close all
end
save([pwd,'/output_data/atom_tracing_model_refinement_total_GPU.mat'],'model_refined')





