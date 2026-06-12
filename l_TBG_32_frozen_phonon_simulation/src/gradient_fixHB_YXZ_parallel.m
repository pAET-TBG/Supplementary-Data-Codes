function [Projs,params,errR,model_arr] = gradient_fixHB_YXZ_parallel(para, xdata, ydata)
fprintf('\nHB gradient algorithm (parallel)\n');

errR = [];
model_arr = [];

Z_arr	   = xdata.Z_arr;
Res        = xdata.Res;
halfWidth  = xdata.halfWidth;
iterations = xdata.iterations;
step_sz    = xdata.step_sz;

model      = xdata.model;
model_ori  = xdata.model_ori;
angles     = xdata.angles;
atom       = xdata.atoms;

num_atom      = numel(atom);
num_atom_type = numel(unique(atom));

[N1,N2,num_pj] = size(ydata);

fixedfa = reshape(make_fixedfa_man([N1 N2], Res, Z_arr), [N1,N2]);

model     = model/Res;
model_ori = model_ori/Res;

[X_crop,Y_crop] = ndgrid(-halfWidth:halfWidth, -halfWidth:halfWidth);
Z_crop = -halfWidth:halfWidth;

h = zeros(1,1,num_atom,'single');
b = zeros(1,1,num_atom,'single');
if size(para,2)==num_atom_type
    for k=1:num_atom_type
        h(atom==k) = para(1,k);
        b(atom==k) = para(2,k);
    end
elseif size(para,2)==num_atom
    h(:) = para(1,:);
    b(:) = para(2,:);
else
    error('para size mismatch');
end
b = (Res*pi)^2./b;

h = single(h);
b = single(b);
X_crop = single(X_crop);
Y_crop = single(Y_crop);
Z_crop = single(Z_crop);

X_ori = reshape(model_ori(1,:), [1,1,num_atom]);
Y_ori = reshape(model_ori(2,:), [1,1,num_atom]);
Z_ori = reshape(model_ori(3,:), [1,1,num_atom]);

X = reshape(model(1,:), [1,1,num_atom]);
Y = reshape(model(2,:), [1,1,num_atom]);
Z = reshape(model(3,:), [1,1,num_atom]);

scale = 1/Res;
N_s = 2*halfWidth+1;

% --- optional: ensure parpool exists (won't error if toolbox missing)
usePar = license('test','Distrib_Computing_Toolbox') && ~isempty(ver('parallel'));
if usePar
    try
        p = gcp('nocreate');
        if isempty(p), parpool; end
    catch
        usePar = false;
    end
end

% Use Constants to reduce broadcast overhead in parfor
if usePar
    c_fixedfa = parallel.pool.Constant(fixedfa);
    c_Xcrop   = parallel.pool.Constant(X_crop);
    c_Ycrop   = parallel.pool.Constant(Y_crop);
    c_Zcrop   = parallel.pool.Constant(Z_crop);
else
    c_fixedfa = [];
    c_Xcrop = []; c_Ycrop = []; c_Zcrop = [];
end

for iter = 1:iterations

    % ---------------------------
    % 1) Forward projection in parallel over i
    % ---------------------------
    Projs_raw = zeros(N1,N2,num_pj,'single');

    % store per-projection per-atom info needed for gradients
    idx_cell = cell(num_pj,1);      % each: int32(2,num_atom)
    gx_cell  = cell(num_pj,1);      % each: single(N_s,N_s,num_atom)
    gy_cell  = cell(num_pj,1);
    gz_cell  = cell(num_pj,1);
    gh_cell  = cell(num_pj,1);
    gb_cell  = cell(num_pj,1);

    if usePar
        parfor i = 1:num_pj
            fixedfa_i = c_fixedfa.Value;
            Xc = c_Xcrop.Value; Yc = c_Ycrop.Value; Zc = c_Zcrop.Value;

            [Proj_i, idx_i, gx_i, gy_i, gz_i, gh_i, gb_i] = ...
                local_forward_oneproj(i, angles, X, Y, Z, h, b, ...
                                      Xc, Yc, Zc, N1, N2, halfWidth, num_atom);

            Projs_raw(:,:,i) = Proj_i;

            idx_cell{i} = idx_i;
            gx_cell{i}  = gx_i;
            gy_cell{i}  = gy_i;
            gz_cell{i}  = gz_i;
            gh_cell{i}  = gh_i;
            gb_cell{i}  = gb_i;
        end
    else
        for i = 1:num_pj
            [Proj_i, idx_i, gx_i, gy_i, gz_i, gh_i, gb_i] = ...
                local_forward_oneproj(i, angles, X, Y, Z, h, b, ...
                                      X_crop, Y_crop, Z_crop, N1, N2, halfWidth, num_atom);

            Projs_raw(:,:,i) = Proj_i;

            idx_cell{i} = idx_i;
            gx_cell{i}  = gx_i;
            gy_cell{i}  = gy_i;
            gz_cell{i}  = gz_i;
            gh_cell{i}  = gh_i;
            gb_cell{i}  = gb_i;
        end
    end

    % ---------------------------
    % 2) Apply fixedfa (can be parallel too)
    % ---------------------------
    Projs = zeros(N1,N2,num_pj,'single');
    if usePar
        parfor i = 1:num_pj
            Projs(:,:,i) = real(my_ifft(my_fft(Projs_raw(:,:,i)) .* c_fixedfa.Value));
        end
    else
        for i = 1:num_pj
            Projs(:,:,i) = real(my_ifft(my_fft(Projs_raw(:,:,i)) .* fixedfa));
        end
    end

    % residual + error
    res = Projs - ydata;
    errR(iter) = sum(abs(Projs(:)-ydata(:))) / sum(abs(ydata(:)));
    fprintf('%d.f = %.5f\n', iter, errR(iter));

    % ---------------------------
    % 3) Gradient accumulation & update X/Y/Z (parallel over atoms k)
    % ---------------------------
    dt = step_sz / mean(h)^2 / mean(b)^2 / halfWidth^2 / num_pj / (N1*N2);

    X_new = X; Y_new = Y; Z_new = Z;

    if usePar
        parfor k = 1:num_atom
            grad_x_k = 0; grad_y_k = 0; grad_z_k = 0;
            % grad_h_k = 0; grad_b_k = 0; % (kept if you want update h,b)

            for i = 1:num_pj
                idx_i = idx_cell{i};
                xi = double(idx_i(1,k));
                yi = double(idx_i(2,k));
                indx = xi + (-halfWidth:halfWidth) + round((N1+1)/2);
                indy = yi + (-halfWidth:halfWidth) + round((N2+1)/2);

                % NOTE: assume no out-of-bound; if needed, add boundary checks here
                rpatch = res(indx, indy, i);

                grad_x_k = grad_x_k + sum(sum(rpatch .* gx_cell{i}(:,:,k)));
                grad_y_k = grad_y_k + sum(sum(rpatch .* gy_cell{i}(:,:,k)));
                grad_z_k = grad_z_k + sum(sum(rpatch .* gz_cell{i}(:,:,k)));

                % grad_h_k = grad_h_k + sum(sum(rpatch .* gh_cell{i}(:,:,k)));
                % grad_b_k = grad_b_k + sum(sum(rpatch .* gb_cell{i}(:,:,k)));
            end

            X_new(k) = X_new(k) - dt * grad_x_k;
            Y_new(k) = Y_new(k) - dt * grad_y_k;
            Z_new(k) = Z_new(k) - dt * grad_z_k;

            % If you want to update h,b, uncomment and tune scaling:
            % b_new(k) = b_new(k) + ...;
            % h_new(k) = h_new(k) - ...;
        end
    else
        for k = 1:num_atom
            grad_x_k = 0; grad_y_k = 0; grad_z_k = 0;

            for i = 1:num_pj
                idx_i = idx_cell{i};
                xi = double(idx_i(1,k));
                yi = double(idx_i(2,k));
                indx = xi + (-halfWidth:halfWidth) + round((N1+1)/2);
                indy = yi + (-halfWidth:halfWidth) + round((N2+1)/2);

                rpatch = res(indx, indy, i);

                grad_x_k = grad_x_k + sum(sum(rpatch .* gx_cell{i}(:,:,k)));
                grad_y_k = grad_y_k + sum(sum(rpatch .* gy_cell{i}(:,:,k)));
                grad_z_k = grad_z_k + sum(sum(rpatch .* gz_cell{i}(:,:,k)));
            end

            X_new(k) = X_new(k) - dt * grad_x_k;
            Y_new(k) = Y_new(k) - dt * grad_y_k;
            Z_new(k) = Z_new(k) - dt * grad_z_k;
        end
    end

    X = X_new; Y = Y_new; Z = Z_new;

    % clamp displacement relative to original (your original logic)
    diff_X = X - X_ori;
    diff_Y = Y - Y_ori;
    diff_Z = Z - Z_ori;
    diff_norm = sqrt(diff_X.^2 + diff_Y.^2 + diff_Z.^2);
    index_3 = diff_norm > scale;

    X(index_3) = X_ori(index_3) + scale * diff_X(index_3) ./ diff_norm(index_3);
    Y(index_3) = Y_ori(index_3) + scale * diff_Y(index_3) ./ diff_norm(index_3);
    Z(index_3) = Z_ori(index_3) + scale * diff_Z(index_3) ./ diff_norm(index_3);

    h = max(h,0);
    b = max(b,0);

    model_arr(:,:,iter) = [X(:),Y(:),Z(:)]' * Res;
end

model  = [X(:),Y(:),Z(:)]' * Res;
params = [h(:)'; (Res*pi)^2 ./ b(:)'; model];

end

% =========================================================================
% Local: compute one projection i (no side effects, safe for parfor)
% =========================================================================
function [Proj_i, idx_i, gx_i, gy_i, gz_i, gh_i, gb_i] = local_forward_oneproj( ...
    i, angles, X, Y, Z, h, b, X_crop, Y_crop, Z_crop, N1, N2, halfWidth, num_atom)

Proj_i = zeros(N1,N2,'single');

RM1 = MatrixQuaternionRot([0 1 0], angles(i,1));
RM2 = MatrixQuaternionRot([1 0 0], angles(i,2));
RM3 = MatrixQuaternionRot([0 0 1], angles(i,3));
R   = RM1*RM2*RM3;

model_rot = R' * [X(:)'; Y(:)'; Z(:)'];

X_cen = reshape(model_rot(1,:), [1,1,num_atom]);
Y_cen = reshape(model_rot(2,:), [1,1,num_atom]);
Z_cen = reshape(model_rot(3,:), [1,1,num_atom]);

X_round = round(X_cen);
Y_round = round(Y_cen);
Z_round = round(Z_cen);

Dx = bsxfun(@plus, X_crop, X_round - X_cen);
Dy = bsxfun(@plus, Y_crop, Y_round - Y_cen);
Dz = bsxfun(@plus, Z_crop, Z_round - Z_cen);

l2_xy = Dx.^2 + Dy.^2;
l2_z  = Dz.^2;

l2_xy_b     = bsxfun(@times, l2_xy, b);
l2_z_b      = bsxfun(@times, l2_z , b);

exp_l2_z_b  = exp(-l2_z_b);
exp_l2_xy_b = exp(-l2_xy_b);

pj_j   = bsxfun(@times, exp_l2_xy_b, sum(exp_l2_z_b));
pj_j   = bsxfun(@times, pj_j, h);
pj_j_b = bsxfun(@times, pj_j, b);

grad_exp = bsxfun(@times, exp_l2_xy_b, sum(l2_z .* exp(-l2_z_b)));
bj_j     = bsxfun(@times, h, grad_exp) + pj_j .* l2_xy;

R2_Dx = (R(1,1)*Dx + R(1,2)*Dy) .* pj_j_b;
R2_Dy = (R(2,1)*Dx + R(2,2)*Dy) .* pj_j_b;
R2_Dz = (R(3,1)*Dx + R(3,2)*Dy) .* pj_j_b;

Dz_exp_l2_z_b = Dz .* exp_l2_z_b;
sum_Dz_exp    = bsxfun(@times, sum(Dz_exp_l2_z_b), exp_l2_xy_b);
sum_Dz_hb     = bsxfun(@times, sum_Dz_exp, h .* b);

xj_j = R2_Dx + R(1,3) * sum_Dz_hb;
yj_j = R2_Dy + R(2,3) * sum_Dz_hb;
zj_j = R2_Dz + R(3,3) * sum_Dz_hb;

idx_i = zeros(2, num_atom, 'int32');

gx_i = zeros(2*halfWidth+1, 2*halfWidth+1, num_atom, 'single');
gy_i = zeros(2*halfWidth+1, 2*halfWidth+1, num_atom, 'single');
gz_i = zeros(2*halfWidth+1, 2*halfWidth+1, num_atom, 'single');
gh_i = zeros(2*halfWidth+1, 2*halfWidth+1, num_atom, 'single');
gb_i = zeros(2*halfWidth+1, 2*halfWidth+1, num_atom, 'single');

for k = 1:num_atom
    xi = X_round(k);
    yi = Y_round(k);

    idx_i(:,k) = int32([xi; yi]);

    indx = xi + (-halfWidth:halfWidth) + round((N1+1)/2);
    indy = yi + (-halfWidth:halfWidth) + round((N2+1)/2);

    Proj_i(indx, indy) = Proj_i(indx, indy) + pj_j(:,:,k);

    gx_i(:,:,k) = xj_j(:,:,k);
    gy_i(:,:,k) = yj_j(:,:,k);
    gz_i(:,:,k) = zj_j(:,:,k);
    gh_i(:,:,k) = pj_j(:,:,k);
    gb_i(:,:,k) = bj_j(:,:,k);
end

end