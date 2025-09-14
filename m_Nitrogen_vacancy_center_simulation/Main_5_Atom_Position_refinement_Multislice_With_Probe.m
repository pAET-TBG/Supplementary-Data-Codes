%% Main Position refinement 
%% care about the edge
% refinement of the atomic coordinates with existed model, type and
% reconstruction volume: minimize the error between the atomic coordinates 
% and the measured projections
clear;clc;

addpath([pwd,'/src/'])
addpath([pwd,'/input_data/'])
addpath([pwd,'/output_data/'])

% add the path to load measured projections and angles, you can comment it 
% and move the projections and angles into input folder
% addpath('../1_Measured_data/') 

% read in files: measured projections and angles; 
% atomic position and types after classification (we use 17 previous)
projections = (importdata('Projs_Tracing_NVC_Multislice_With_Probe_Total.mat'));
angles      = (importdata('angles_Tracing_NVC_Multislice_With_Probe_Total.mat'));
model       = (importdata('atom_Tracing_NVC_Multislice_With_Probe_Total.mat'));
model_org   = (importdata('atom_Tracing_NVC_Multislice_With_Probe_Total.mat'));
atoms       = (importdata('label_Tracing_NVC_Multislice_With_Probe_Total.mat'));

% process the data: has to be double type in lsqcurvefit;
%% the projections dimension has to be odd
projections = max(projections,0);                                          % the positivity constrain for each projections
projections = projections(1:end,1:end,:);
projections = My_paddzero(projections,size(projections) + [50 50 0],'double');

[N1,N2,num_pj] = size(projections);
% the cropped bondary size for each atoms
halfWidth = 10;
% the atomic number for different type:
% use 28 for type 1, 45 for type 2, 78 for type 3; (C: 6 N: 7)
Z_arr   = 6;     
% indicate the pixel size for measured projections
Res     = 0.2129;            

xdata = [];
xdata.Res       = Res;                                                     % the resolution of the pixel for each projections
xdata.Z_arr     = Z_arr;                                                   % the values for different type of atoms 
xdata.halfWidth = halfWidth;                                               % the cropped bondary size for each atoms
xdata.atoms     = atoms;                                                   % the types of each atoms
xdata.model     = model;                                                   % the positions of each atoms
xdata.angles    = angles;                                                  % the tilting angle for each projections
para0 = [1 ;  
         4.5];                                                             % the initial parameter for the fitted function
lb = [0.1;
      1];                                                                  % the lower bound of each paremeter
ub = [2 ; 
      7];                                                                  % the upper bound of each parameter
model_refined = model;                                                     % the positions of each atom prepare for refining

% option method for optimization
opt = optimset('TolFun', 1e-12, 'TolX', 1e-8, 'MaxIter', 1000, 'Display', 'iter');

position_grad = zeros(size(model));

for jjjj=1:10
    fprintf('iteration num: %d; \n',jjjj);
    x0 = para0;                                                            % the parameter of the function perpared for optimization
    x0(1,:)=x0(1,:)/x0(1,1);                                               % normalize the para0(1,1) para0(1,2) para0(1,3)
    
    xdata.model = model;                                                   % the positions of each atoms perpared for optimization
    xdata.model_ori = model;                                               % the original position of each atoms will not be changed during this iteration
    xdata.projections=projections;
    
    % Cal_Bproj_2tpue2 function is used to calculate the projecitons according to the atom positions model
    % @fun = @Cal_Bproj_2type2 is the function to calculate the projections
    % x0 = x0 is the parameter for Cal_Bproj_2tpye2
    % xdata = xdata is the input value for the function Cal_Bproj_2type2
    % ydata = projections is the value Cal_Bproj_2type2 needs to match
    % lb = lb is the low boundary for x0
    % ub = ub is the up boundary for x0
    % opt = opt is the option method for optimization
    % to obtain the most match p7arameter in this iteration
    
    [para0, ~,~] = lsqcurvefit(@Cal_Bproj_2type2, x0, xdata, projections, lb, ub, opt);
    [y_pred,~] = Cal_Bproj_2type(para0, xdata, projections);
    
    xdata.projections = [];
    xdata.step_sz    = 1;
    xdata.iterations = 10;
    [y_pred,para0,errR] = gradient_B_2type_difB(para0, xdata, projections);
    
    xdata.step_sz    = 1;
    xdata.iterations = 10;
    [y_pred,para,errR] = gradient_fixHB_XYZ(para0, xdata, projections); 
    model_refined = para(3:5,:);
    
    position_grad = 1.*(model_refined - model) + 0.0.* position_grad;
    model = model + position_grad;
end
save([pwd,'/output_data/atom_Tracing_NVC_Multislice_With_Probe_Tracing_Refinement.mat'],'model')

%% calculate the B factor
calculate_B_factor(para0, xdata, model, projections);

