function [Projs,param] = Cal_Bproj_2type2(para, xdata)
% tic
para=abs(para);                                                            % the parameter's positivity constrain (parameters have to be positive)

Z_arr	   = xdata.Z_arr;                                                  % use 28 for type 1, 45 for type 2, 78 for type 3 
                                                                           % the values for different type of atoms (Maybe it is the Sigma(standard deviation) of Gaussian)
Res        = xdata.Res;                                                    % the resolution of the pixel for each projections     
halfWidth  = xdata.halfWidth;                                              % the cropped bondary size for each atoms

model      = xdata.model;                                                  % the positions of each atoms
angles     = xdata.angles;                                                 % the tilting angles for each projections
atom       = xdata.atoms;                                                  % the types of each atoms
num_atom   = numel(atom);                                                  % the number of total atoms
atom_type_num = numel(unique(atom));                                       % the number of atom types
ydata      = xdata.projections;                                            % the ground truth projections 

[N1,N2,~] = size(ydata);                                                   % obtain the number of column (N1); and the number of row (N2)
num_pj=size(angles,1);                                                     % obtain the number of ptojections
% N_s = 2*halfWidth+1;

%%%%%%% (fixefa might be a filter in the reciprocal space) %%%%%%%%%%
fixedfa = reshape( make_fixedfa_man([N1,N2], Res, Z_arr), [N1,N2]); %       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% max_fa = max(abs(fixedfa(:)));
model = model/Res;                                                         % the model's unit from angstrom to pixel size

dtype = 'single';
X_rot = zeros(num_pj,num_atom,dtype);                                      % initialize the matrix for the X (column) positions after rotation
Y_rot = zeros(num_pj,num_atom,dtype);                                      % initialize the matrix for the Y (row) positions after rotation
Z_rot = zeros(num_pj,num_atom,dtype);                                      % initialize the matrix for the Z positions after rotation

% calculate the atom positions after rotation
for i=1:num_pj
  R1 = MatrixQuaternionRot([0 0 1],angles(i,1));  
  R2 = MatrixQuaternionRot([0 1 0],angles(i,2));
  R3 = MatrixQuaternionRot([1 0 0],angles(i,3));  
  R   = (R1*R2*R3)';
      
  rotCoords = R*model;
  X_rot(i,:) = rotCoords(1,:);
  Y_rot(i,:) = rotCoords(2,:);
  Z_rot(i,:) = rotCoords(3,:);
end

% generate the coordination for the crop region in the 3D model
[X_crop,Y_crop] = ndgrid( -halfWidth:halfWidth, -halfWidth:halfWidth);      
Z_crop = -halfWidth:halfWidth;

% normalize the parameter and adjust the parameter
para = reshape(para,[2 atom_type_num]);
h       = para(1,:)/para(1,1);                                             % the height of the Gaussian fitting
b       = (pi*Res)^2 ./ para(2,:);                                         % the width of the Gaussian fitting (Standard Deviation)

% obtain the total number of each type of atoms
num_atom_type = zeros(atom_type_num,1);
for j=1:atom_type_num
    num_atom_type(j) = nnz(atom==j);
end

Grad    = zeros(N1,N2,num_pj, 3);
for i=1:num_pj
    for j=1:atom_type_num
        atom_type_j = atom==j;
        
        X_cen = reshape(X_rot(i,atom_type_j), [1,1,num_atom_type(j)]);
        Y_cen = reshape(Y_rot(i,atom_type_j), [1,1,num_atom_type(j)]);
        Z_cen = reshape(Z_rot(i,atom_type_j), [1,1,num_atom_type(j)]);
        
        X_round = round(X_cen);
        Y_round = round(Y_cen);
        Z_round = round(Z_cen);
        
        l2_xy = bsxfun(@plus, X_crop, X_round - X_cen).^2 + ...
                bsxfun(@plus, Y_crop, Y_round - Y_cen).^2;
        
        l2_z  = bsxfun(@plus, Z_crop, Z_round - Z_cen).^2    ;        
        
        pj_j   = bsxfun(@times, exp(-l2_xy* b(j) ), sum(exp(-l2_z* b(j))) );
        pj_j_h = h(j)*pj_j;
        
        for k=1:num_atom_type(j)
            indx = X_round(k) + (-halfWidth:halfWidth) + round((N1+1)/2);
            indy = Y_round(k) + (-halfWidth:halfWidth) + round((N2+1)/2);

            Grad(indx,indy,i,j)   = Grad(indx,indy,i,j)   + pj_j_h(:,:,k);
        end
    end
end

Projs = sum(Grad,4);

for i=1:num_pj
    Projs(:,:,i) = (( my_ifft( my_fft(Projs(:,:,i)) .* (fixedfa) ) ));
end

k=sum(Projs(:).*ydata(:))/sum(Projs(:).^2);
Projs=Projs*k;
param=[k*h; (pi*Res)^2 ./ b(:)'];
end