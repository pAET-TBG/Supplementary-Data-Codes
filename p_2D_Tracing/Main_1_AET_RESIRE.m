clear
clc
addpath([pwd,'/src/'])
addpath([pwd,'/gpu_src/'])
addpath([pwd,'/input_data/']);

%% load data
projections  = importdata([pwd,'/input_data/Projections.mat']);
angles       = importdata([pwd,'/input_data/Angles.mat' ]);
Mask         = importdata([pwd,'/input_data/Mask.mat' ]);

%% pre-process: rotation setupx
rotation       = 'ZYX';  % Euler angles setting ZYZ
dtype          = 'single';
projections_refined = cast(projections,dtype);
angles_refined      = cast(angles,dtype);
clear projections angles

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
% extract size of projections & num of projections
[dimx, dimy, Num_pj] = size(projections_refined);
dimz = 80;
vec1 = matR(:,1);
vec2 = matR(:,2);
vec3 = matR(:,3);

%% generate rotation matrices
Rs = zeros(3,3,Num_pj, dtype);
for k = 1:Num_pj
    phi   = angles_refined(k,1);
    theta = angles_refined(k,2);
    psi   = angles_refined(k,3);
    
    % compute rotation matrix R w.r.t euler angles {phi,theta,psi}
    rotmat1 = MatrixQuaternionRot(vec1,phi);
    rotmat2 = MatrixQuaternionRot(vec2,theta);
    rotmat3 = MatrixQuaternionRot(vec3,psi);
    R =  single(rotmat1*rotmat2*rotmat3)';
    Rs(:,:,k) = R;
end

%% parameter
step_size      = 0.1;  %step_size <=1 but can be larger is sparse
iterations     = 200;
positivity     = true;

%% Test resire2 code on thinfilm with multi-GPU
fprintf('\nResire code:thinfilm (multi-GPU)\n');
dim_ext = [dimx, dimy, dimz];
[rec] = RT3_film_multiGPU( (projections_refined), (Rs), dim_ext, ...
    (iterations), (step_size) , (positivity) );
rec   = rec.*Mask;
save([pwd,'/output_data/tBLG_reconstruction_volume.mat'],'rec')



