clear
clc
addpath([pwd,'/src/'])
addpath([pwd,'/gpu_src/'])
addpath([pwd,'/cuda_code/'])
addpath([pwd,'/input_data/']);

%% load data
projections  = importdata([pwd,'/input_data/Projections.mat']);
angles       = importdata([pwd,'/input_data/Angles.mat' ]);

%% input
rotation       = 'ZYX';  % Euler angles setting ZYZ
dtype          = 'single';
projections    = cast(projections,dtype);
defocus_param  = cast(zeros(size(angles,1),1),dtype);
angles         = cast(angles,dtype);

% compute normal vector of rotation matrix
matR = zeros(3,3);
if length(rotation)~=3
    disp('rotation not recognized. Set rotation = ZYX\n'); rotation = 'ZYX';
end
for i=1:3
    switch rotation(i)
        case 'X',   matR(:,i) = [1;0;0];
        case 'Y',   matR(:,i) = [0;1;0];
        case 'Z',   matR(:,i) = [0;0;1];
        otherwise,  matR = [0,0,1;
                0,1,0;
                1,0,0];
            disp('Rotation not recognized. Set rotation = ZYX');
            break
    end
end
vec1 = matR(:,1); vec2 = matR(:,2); vec3 = matR(:,3);

% extract size of projections & num of projections
[dimx, dimy, Num_pj] = size(projections);

%% rotation matrix
Rs = zeros(3,3,Num_pj, dtype);
for k = 1:Num_pj
    phi   = angles(k,1);
    theta = angles(k,2);
    psi   = angles(k,3);
    
    % compute rotation matrix R w.r.t euler angles {phi,theta,psi}
    rotmat1 = MatrixQuaternionRot(vec1,phi);
    rotmat2 = MatrixQuaternionRot(vec2,theta);
    rotmat3 = MatrixQuaternionRot(vec3,psi);
    R =  single(rotmat1*rotmat2*rotmat3)';
    Rs(:,:,k) = R;
end

%% parameters
dimz           = 80;
positivity     = 1;                                                        % 1 means true, 0 means false
l2_regularizer = 0.00;                                                     % a small positive number
is_avg_on_y    = 0;                                                        % 1 means averaging in the y-direction, 0 means no

defocus_step   = 80;     
semi_angle     = 32.0e-3;    
Voltage        = 70*10^3;
pixelsize      = 0.19/2;                            
nr             = 200; 
defocus_scale  = 1;

% make sure these following parameters are double
defocus_info = [Voltage, pixelsize, nr, semi_angle, defocus_step, defocus_scale];
constraints  = [positivity, is_avg_on_y, l2_regularizer];

% only projections and Rs are single
%% generate support
% first time reconstruction to figure out the two peak
rec        = zeros([dimx,dimy,dimz]);
step_size  = 1;                                                            %step_size <=1 but can be larger is sparse
iterations = 200; 
GD_info    = [iterations, step_size];
[rec, cal_proj4] = RT3_defocus_1GPU( (projections), (Rs), (dimz), GD_info , (constraints),  ...
        defocus_info, defocus_param,rec);

% crop the reconstruction region without contamination
rec_crop = rec(100:end-100,100:end-100,:);

% set the threshold and average the result of different threshold
precentageVector = 0.40:0.0001:0.60;

sumzVector       = zeros(dimz,size(precentageVector,2));
count            = 1;
for pre = precentageVector
    rec_crop(rec_crop<max(rec_crop(:))*pre) = 0;
    sumz = sum(sum(rec_crop,1),2);
    sumz = sumz - min(sumz);sumz = sumz./max(sumz);
    sumzVector(:,count) = sumz;
    count = count + 1;
end
sumzAve = mean(sumzVector,2);
[minVal,ind] = min(sumzAve(1:19));
sumzAve(1:ind) = minVal;
[minVal,ind] = min(sumzAve(62:end));
sumzAve(61+ind:end) = minVal;
peak(2) = fit_Guassian_function(sumzAve(56:64),56:64,0);
peak(1) = fit_Guassian_function(sumzAve(19:27),19:27,0);

%% generalized Gauss
beta =4;
p0 = [1 ,peak(1) ,10 ,1 ,peak(2), 10];
generalizedGauss = @(p,x) p(1).*exp(-abs(((x-p(2))./(p(3)))).^beta) + p(4).*exp(-abs(((x-p(5))./(p(6)))).^beta);
gaussian_curve = generalizedGauss(p0,1:dimz);
gaussian_curve = gaussian_curve - min(gaussian_curve(:));
gaussian_curve = gaussian_curve./max(gaussian_curve(:)).*0.01;           % for Gaussian it is 0.01
gaussian_curve = gaussian_curve + 1 - max(gaussian_curve(:));
rec = zeros([dimx,dimy,dimz]);
support = zeros(size(rec));
for tempz = 1:dimz
    support(:,:,tempz) = gaussian_curve(tempz);
end

%% iteration: minimize ||Au-b||^2 by gradient descent: run reconstruction the first time
% syntax 1: no initial rec
step_size  = 0.1;  %step_size <=1 but can be larger is sparse
iterations = 1;
GD_info    = [iterations, step_size];
for iter = 1:200
    tic
    rec = rec.*support;
    [rec, cal_proj4] = RT3_defocus_1GPU( (projections), (Rs), (dimz), GD_info , (constraints),  ...
        defocus_info, defocus_param,rec);
    toc
end
%% store the large field of view reconstruction result error =0.0609547
save([pwd,'/output_data/tBLG_reconstruction_volume.mat'],'rec')





