clear
clc
addpath([pwd,'/src/'])
addpath([pwd,'/gpu_src/'])
addpath([pwd,'/cuda_code/'])
addpath([pwd,'/input_data/']);

%% load data
projections = importdata([pwd '/input_data/Projs_NVCenter_Multislice.mat' ]);
angles      = importdata([pwd '/input_data/Angles_NVcenter.mat' ]);

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
dimz           = 800;
positivity     = 1;         % 1 means true, 0 means false
l2_regularizer = 0.00;      % a small positive number
is_avg_on_y    = 0;         % 1 means averaging in the y-direction, 0 means no

defocus_step   = 800;     
semi_angle     = 32.0e-3;    
Voltage        = 70*10^3;
pixelsize      = 0.2129;                      
nr             = 200; 
defocus_scale  = 1;

% make sure these following parameters are double
defocus_info = [Voltage, pixelsize, nr, semi_angle, defocus_step, defocus_scale];
constraints  = [positivity, is_avg_on_y, l2_regularizer];

% only projections and Rs are single
rec        = zeros([dimx,dimy,dimz]);
step_size  = 1;  %step_size <=1 but can be larger is sparse
iterations = 1000; 
GD_info    = [iterations, step_size];
[rec, cal_proj4] = RT3_defocus_1GPU( (projections), (Rs), (dimz), GD_info , (constraints),  ...
        defocus_info, defocus_param,rec);

%% store the large field of view reconstruction result error = 0.0202654
save([pwd,'/output_data/Rec_With_Probe_NVcenter_Multislice.mat'],'rec','-v7.3')
