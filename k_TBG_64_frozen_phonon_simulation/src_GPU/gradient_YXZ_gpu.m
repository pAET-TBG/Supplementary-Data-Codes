function [params,errR,model_out, proj_iter,grad_x_set,grad_y_set,grad_z_set] =...
    gradient_YXZ_gpu(para, xdata, mea_proj)
fprintf('\nHB gradient algorithm\n');
errR=[];
Z_arr	   = xdata.Z_arr;
Res        = xdata.Res;
halfWidth  = xdata.halfWidth;
iterations = xdata.iterations;
step_sz    = xdata.step_sz;

model      = xdata.model;
angles     = xdata.angles;
atoms       = xdata.atoms;
num_atom   = numel(atoms);
num_atom_type = numel(unique(atoms));

[N1,N2,num_pj] = size(mea_proj);

% atom_type1 = xdata.atom_type1;
% atom_type2 = xdata.atom_type2;
% atom_type3 = xdata.atom_type3;
% X_rot = xdata.X_rot;
% Y_rot = xdata.Y_rot;
% Z_rot = xdata.Z_rot;

fixedfa = reshape( make_fixedfa_man([N1 N2], Res, Z_arr), [N1,N2] );
model = model/Res;
% model_ori = model_ori/Res;
%
% for i=1:num_pj
%     ydata(:,:,i) = max(real( my_ifft( my_fft( ydata(:,:,i) ) ./fixedfa )),0);
% end

% l2_type1 =  xdata.l2_type1;
% l2_type2 =  xdata.l2_type2;
% l2_type3 =  xdata.l2_type3;

% h1 = para(1,1);
% h2 = para(1,2);
% h3 = para(1,3);
%
% b1 = para(2,1);
% b2 = para(2,2);
% b3 = para(2,3);


% [X_crop,Y_crop] = ndgrid( -halfWidth:halfWidth, -halfWidth:halfWidth);
% Z_crop = -halfWidth:halfWidth;

%Y_crop = reshape(Y_crop, [1,2*halfWidth+1,2*halfWidth+1]);
%X_crop = reshape(X_crop, [1,2*halfWidth+1,2*halfWidth+1]);
%Z_crop = reshape(Z_crop, [1,2*halfWidth+1]);

h_all = zeros(1,1,num_atom, 'single');
b_all = zeros(1,1,num_atom, 'single');
if size(para,2)==num_atom_type
    for k=1:num_atom_type
        h_all(atoms==k) = para(1,k);
        b_all(atoms==k) = para(2,k);
    end
elseif size(para,2)==num_atom
    h_all(:) = para(1,:);
    b_all(:) = para(2,:);
else
    display('error')
    return
end
b_all = (Res*pi)^2./b_all;

h       = para(1,:);
b       = (pi*Res)^2 ./ para(2,:);

h = (single(h));
b = (single(b));



% N_s = 2*halfWidth+1;

% grad_h_set = zeros(N_s,N_s,num_pj,num_atom,'single');
% grad_b_set = zeros(N_s,N_s,num_pj,num_atom,'single');
% grad_x_set = zeros(N_s,N_s,num_pj, num_atom,'single');
% grad_y_set = zeros(N_s,N_s,num_pj, num_atom,'single');
% grad_z_set = zeros(N_s,N_s,num_pj, num_atom,'single');


% X_crop = (single(X_crop));
% Y_crop = (single(Y_crop));
% Z_crop = (single(Z_crop));
% 
% X_ori = reshape( model_ori(1,:), [1,1,num_atom] );
% Y_ori = reshape( model_ori(2,:), [1,1,num_atom] );
% Z_ori = reshape( model_ori(3,:), [1,1,num_atom] );
% X     = reshape( model(1,:), [1,1,num_atom] );
% Y     = reshape( model(2,:), [1,1,num_atom] );
% Z     = reshape( model(3,:), [1,1,num_atom] );
scale=1/Res;
dt = step_sz/mean(h_all)^2/mean(b_all)^2/halfWidth^2/num_pj/(N1*N2);
fprintf("dt=%.12f,%f,%f,%f\n", dt, step_sz, mean(h), mean(b));

Rs = zeros(3,3,num_pj,'single');
for i=1:num_pj
    RM1 = MatrixQuaternionRot([0 1 0],angles(i,1));
    RM2 = MatrixQuaternionRot([1 0 0],angles(i,2));
    RM3 = MatrixQuaternionRot([0 0 1],angles(i,3));
    Rs(:,:,i) = single(RM1*RM2*RM3)';
end
disp(model(:,1))

    % model    = gpuArray(single(rand(3, num_atoms)));       % [3, num_atoms]
    % atoms    = gpuArray(int32(randi([1 num_types], 1, num_atoms))); % [1, num_atoms]
    % mea_proj = gpuArray(single(rand(N1, N2, num_pj)));     % [N1, N2, num_pj]
    % b        = gpuArray(single(rand(1, num_types)));       % [num_types]
    % h        = gpuArray(single(rand(1, num_types)));       % [num_types]
    % Rs       = gpuArray(single(rand(3, 3, num_pj)));       % [3,3,num_pj]

    % Call the MEX function (assuming compiled as accumulateProjection.mex*)
    % Grad_model
[model_iter, proj_iter, grad_x_set, grad_y_set, grad_z_set] = ...
    cal_grad_all(single(model), int32(atoms), single(mea_proj), single(fixedfa), ...
    single(b), single(h), single(Rs), single(scale), single(iterations), single(dt), halfWidth);

    % model_arr(:,:,iter) = model_iter*Res;
%%
% for i=1:num_pj
%     Projs(:,:,i) = max( 0, real( my_ifft( my_fft(Projs(:,:,i)) .* fixedfa ) ) );
% end
model_out = model_iter*Res;
params = [h(:)'; (Res*pi)^2./b(:)'];
% grad_x_set = permute(grad_x_set,[1,2,4,3]);
% grad_y_set = permute(grad_y_set,[1,2,4,3]);
% grad_z_set = permute(grad_z_set,[1,2,4,3]);
end