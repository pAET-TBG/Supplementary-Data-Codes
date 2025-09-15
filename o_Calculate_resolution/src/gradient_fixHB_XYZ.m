function [Projs,params,errR,model_arr] = gradient_fixHB_XYZ(para, xdata, ydata)
fprintf('\nHB gradient algorithm\n');
errR=[];
model_arr=[];
Z_arr	   = xdata.Z_arr;
Res        = xdata.Res;
halfWidth  = xdata.halfWidth;
iterations = xdata.iterations;
step_sz    = xdata.step_sz;

model      = xdata.model;
model_ori  = xdata.model_ori;
angles     = xdata.angles;
atom       = xdata.atoms;
num_atom   = numel(atom);
num_atom_type = numel(unique(atom));

[N1,N2,num_pj] = size(ydata);

fixedfa = reshape( make_fixedfa_man([N1 N2], Res, Z_arr), [N1,N2] );
model = model/Res;
model_ori = model_ori/Res;

[X_crop,Y_crop] = ndgrid( -halfWidth:halfWidth, -halfWidth:halfWidth);
Z_crop = -halfWidth:halfWidth;

h = zeros(1,1,num_atom, 'single');
b = zeros(1,1,num_atom, 'single');
if size(para,2)==num_atom_type
    for k=1:num_atom_type
        h(atom==k) = para(1,k);
        b(atom==k) = para(2,k);
    end
elseif size(para,2)==num_atom
    h(:) = para(1,:);
    b(:) = para(2,:);
else
    display('error')
    return
end
b = (Res*pi)^2./b;

N_s = 2*halfWidth+1;
index = zeros(2,num_pj, num_atom);
grad_h_set = zeros(N_s,N_s,num_pj,num_atom,'single');
grad_b_set = zeros(N_s,N_s,num_pj,num_atom,'single');
grad_x_set = zeros(N_s,N_s,num_pj, num_atom,'single');
grad_y_set = zeros(N_s,N_s,num_pj, num_atom,'single');
grad_z_set = zeros(N_s,N_s,num_pj, num_atom,'single');

h = (single(h));
b = (single(b));
X_crop = (single(X_crop));
Y_crop = (single(Y_crop));
Z_crop = (single(Z_crop));

X_ori = reshape( model_ori(1,:), [1,1,num_atom] );
Y_ori = reshape( model_ori(2,:), [1,1,num_atom] );
Z_ori = reshape( model_ori(3,:), [1,1,num_atom] );
X     = reshape( model(1,:), [1,1,num_atom] );
Y     = reshape( model(2,:), [1,1,num_atom] );
Z     = reshape( model(3,:), [1,1,num_atom] );
scale=1/Res;

for iter=1:iterations

    Projs = zeros(N1,N2,num_pj,'single');
    for i=1:num_pj
        
        RM1 = MatrixQuaternionRot([0 0 1],angles(i,1));
        RM2 = MatrixQuaternionRot([0 1 0],angles(i,2));
        RM3 = MatrixQuaternionRot([1 0 0],angles(i,3));
        R   = RM1*RM2*RM3;
        model_rot = R'*[X(:)';Y(:)';Z(:)'];
        
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
        
        pj_j     = bsxfun(@times, exp_l2_xy_b, sum(exp_l2_z_b) );
        pj_j     = bsxfun(@times, pj_j , h);
        pj_j_b   = bsxfun(@times, pj_j , b);
        
        grad_exp = bsxfun(@times, exp_l2_xy_b, sum(l2_z.*exp(-l2_z_b)) );
        bj_j     = bsxfun(@times, h, grad_exp) + pj_j.*l2_xy ;
        
        R2_Dx = ( R(1,1)*Dx + R(1,2)*Dy ) .* pj_j_b;
        R2_Dy = ( R(2,1)*Dx + R(2,2)*Dy ) .* pj_j_b;
        R2_Dz = ( R(3,1)*Dx + R(3,2)*Dy ) .* pj_j_b;
        
        Dz_exp_l2_z_b = Dz .* exp_l2_z_b;
        sum_Dz_exp    = bsxfun(@times, sum( Dz_exp_l2_z_b ), exp_l2_xy_b);
        
        sum_Dz_hb     = bsxfun(@times, sum_Dz_exp, h.*b);        
        
        xj_j = R2_Dx + R(1,3) * sum_Dz_hb;
        yj_j = R2_Dy + R(2,3) * sum_Dz_hb;
        zj_j = R2_Dz + R(3,3) * sum_Dz_hb;
        
        for k=1:num_atom
            indx = X_round(k) + (-halfWidth:halfWidth) + round((N1+1)/2);
            indy = Y_round(k) + (-halfWidth:halfWidth) + round((N2+1)/2);
            Projs(indx,indy,i) = Projs(indx,indy,i) + pj_j(:,:,k);
            index(:,i,k) = [X_round(k);Y_round(k)];
            
            grad_h_set(:,:,i,k) = pj_j(:,:,k);
            grad_b_set(:,:,i,k) = bj_j(:,:,k);
            
            grad_x_set(:,:,i,k) = xj_j(:,:,k);
            grad_y_set(:,:,i,k) = yj_j(:,:,k);
            grad_z_set(:,:,i,k) = zj_j(:,:,k);
        end
    end
    
    
    for i=1:num_pj
        Projs(:,:,i) = real(( my_ifft( my_fft(Projs(:,:,i)) .* fixedfa ) ));
    end
    res = Projs - ydata;
    errR(iter)=sum( abs(Projs(:)-ydata(:)) ) / sum(abs(ydata(:)));
    
    for k=1:num_atom
        grad_x_k=0; grad_y_k=0; grad_z_k=0;
        grad_h_k=0; grad_b_k=0;
        for i=1:num_pj
            indx = index(1,i,k) + (-halfWidth:halfWidth) + round((N1+1)/2);
            indy = index(2,i,k) + (-halfWidth:halfWidth) + round((N2+1)/2);
            grad_x_k = grad_x_k + sum(sum( res(indx,indy,i) .* grad_x_set(:,:,i,k) ));
            grad_y_k = grad_y_k + sum(sum( res(indx,indy,i) .* grad_y_set(:,:,i,k) ));
            grad_z_k = grad_z_k + sum(sum( res(indx,indy,i) .* grad_z_set(:,:,i,k) ));
            
            grad_h_k = grad_h_k + sum(sum( res(indx,indy,i) .* grad_h_set(:,:,i,k) ));
            grad_b_k = grad_b_k + sum(sum( res(indx,indy,i) .* grad_b_set(:,:,i,k) ));
        end
        dt = step_sz/mean(h)^2/mean(b)^2/halfWidth^2/num_pj/(N1*N2);
        
        X(k) = X(k) - dt*grad_x_k;
        Y(k) = Y(k) - dt*grad_y_k;
        Z(k) = Z(k) - dt*grad_z_k;
        
        % choose wheter to update h and b
%         b(k) = b(k) + t/halfWidth^4/mean(h)^2 * grad_b_k ;
%         h(k) = h(k) - t/mean(h) * grad_h_k;
    end
    diff_X = X-X_ori;
    diff_Y = Y-Y_ori;
    diff_Z = Z-Z_ori;
    diff_norm = sqrt( diff_X.^2 + diff_Y.^2 + diff_Z.^2 );
    index_3 = diff_norm>scale;
    X(index_3) = X_ori(index_3) + scale*diff_X(index_3) ./ diff_norm(index_3);
    Y(index_3) = Y_ori(index_3) + scale*diff_Y(index_3) ./ diff_norm(index_3);
    Z(index_3) = Z_ori(index_3) + scale*diff_Z(index_3) ./ diff_norm(index_3);
     
    h = max(h,0);
    b = max(b,0);
    
    model_arr(:,:,iter) = [X(:),Y(:),Z(:)]'*Res;
end
model = [X(:),Y(:),Z(:)]'*Res;
params = [h(:)'; (Res*pi)^2./b(:)'; model];
end