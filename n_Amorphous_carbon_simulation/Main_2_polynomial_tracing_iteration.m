%% Main Polynomial Tracing
% code to perform polynomial tracing on reconstruction volume
% reconstruction volume should be upsampled by 3*3*3 linear interpolation
% each local maximum is fitted with 9*9*9 voxels (3*3*3 before interpolation)
% by 4th order polynomial equation to get the position

clear;clc;
addpath('src/');
% add the path to load reconstruction volume, you can comment it and move
% the reconstruction into input folder
% load reconstructed 3D volume
Dsetvol0 = importdata([pwd,'/output_data/tBLG_reconstruction_volume.mat']);
Mask     = importdata([pwd,'/input_data/Mask.mat']);

Dsetvol=Mask.*Dsetvol0;

N = size(Dsetvol);
Dsetvol0 = My_stripzero(Dsetvol,[N(1),N(2),N(3)]);
Dsetvol_1 = Dsetvol0(:,1:end,:);

FinalVol_tracing = (Dsetvol_1);
N1=size(FinalVol_tracing);

%% oversamplling the reconstruction matrix by 3*3*3 by linear interpolation
ov=3;% oversamplling
xx = (1:size(FinalVol_tracing,1)) - round((size(FinalVol_tracing,1)+1)/2);
yy = (1:size(FinalVol_tracing,2)) - round((size(FinalVol_tracing,2)+1)/2);
zz = (1:size(FinalVol_tracing,3)) - round((size(FinalVol_tracing,3)+1)/2);

xxi = ((ov*xx(1)):(xx(end)*ov))/ov;
yyi = ((ov*yy(1)):(yy(end)*ov))/ov;
zzi = ((ov*zz(1)):(zz(end)*ov))/ov;

xxi = xxi(3:end); yyi = yyi(3:end); zzi = zzi(3:end);

[Y,X,Z]     = meshgrid(yy,xx,zz);
[Yi,Xi,Zi]  = meshgrid(yyi,xxi,zzi);

FinalVol_tracing=single(FinalVol_tracing);
FinalVol_tracing_interp     = interp3(Y,X,Z,FinalVol_tracing,Yi,Xi,Zi,'spline',0);
FinalVol    = My_paddzero(FinalVol_tracing_interp,size(FinalVol_tracing_interp)+20);


FinalVol=double(FinalVol);
N2=size(FinalVol);

%%  get polynomial power array
[ii,jj,kk] = meshgrid(0:4,0:4,0:4);
fitCoeff = [ii(:),jj(:),kk(:),0*kk(:)];
fitCoeff(sum(fitCoeff,2)>4,:) = [];
fitCoeff(max(fitCoeff,[],2) == 4,4) = -1;
% return
%%
parpool_size=8;
if parpool_size~=0
    pjob = gcp('nocreate');
    if isempty(pjob)
        parpool(parpool_size)
    elseif pjob.NumWorkers ~= parpool_size
        delete(pjob)
        parpool(parpool_size)
    end
end
%% Tracing parameters

MaxIter = 14;          CritIter = 7;           Th        = 0.00;  
Res     = 0.33967/ov;  minDist  = 1.0/Res;     SearchRad = 3;% 

%% get the local maxima from the reconstruction volume    
se           = strel3d(3);
dilatedBW    = imdilate(FinalVol, se);
maxPos       = find( FinalVol == dilatedBW & FinalVol > Th );
maxVals      = FinalVol(maxPos);
[~,sortInd]  = sort(maxVals,'descend');

maxNum      = min(100000,length(sortInd));
maxPos       = maxPos(sortInd(1:maxNum));
    
fprintf(1,'numpeak = %d \n',length(maxPos));
 
maxXYZ = zeros(length(maxPos),3);
for i=1:length(maxPos)
    [xx,yy,zz] = ind2sub(size(FinalVol),maxPos(i));
    maxXYZ(i,:) = [xx yy zz];
end

clear Xi Yi Zi Dsetvol dilatedBW
    %%
Q = 0.5;  Alpha = 1;
cropHalfSize = SearchRad;
[X,Y,Z] = ndgrid(-cropHalfSize:cropHalfSize,-cropHalfSize:cropHalfSize,-cropHalfSize:cropHalfSize);

SphereInd = find( X.^2+Y.^2+Z.^2 <= (SearchRad+0.5)^2 );
X_sph = X(SphereInd);
Y_sph = Y(SphereInd);
Z_sph = Z(SphereInd);
    
XYZdata.X = X_sph;
XYZdata.Y = Y_sph;
XYZdata.Z = Z_sph;
  
Orders = fitCoeff(:,1:3);
PosArr = zeros(size(maxXYZ));
TotPosArr = zeros(size(maxXYZ)); 
exitFlagArr = zeros(1, size(maxXYZ,1)); 
CoeffArr = repmat(fitCoeff(:,4),[1 size(maxXYZ,1)]);

parfor i=1:size(maxXYZ,1)
    tic
    endFlag = 0;
    consecAccum = 0;
    iterNum = 0;
    while ~endFlag
        iterNum = iterNum + 1;
        if iterNum>MaxIter
            exitFlagArr(i) = -4;
            endFlag = 1;
        end
        
        cropXind = maxXYZ(i,1) + (-cropHalfSize:cropHalfSize);
        cropYind = maxXYZ(i,2) + (-cropHalfSize:cropHalfSize);
        cropZind = maxXYZ(i,3) + (-cropHalfSize:cropHalfSize);
        
        cropVol  = FinalVol(cropXind,cropYind,cropZind);
        
        Pos = PosArr(i,:);
        GaussWeight = exp(-1*Alpha*( (X_sph-Pos(1)).^2 + (Y_sph-Pos(2)).^2 + (Z_sph-Pos(3)).^2 ) / cropHalfSize^2 ); 
        fun = @(p,xdata) calculate_3D_polynomial_Rogers(xdata.X,xdata.Y,xdata.Z,Pos,Orders,p).*GaussWeight; 
        opts = optimset('Display','off');
        [p1,fminres1] = lsqcurvefit(fun,CoeffArr(:,i),XYZdata,cropVol(SphereInd).*GaussWeight,[],[],opts);
        CoeffArr(:,i) = p1;
        %%
        [dX, dY, dZ] = calc_dX_dY_dZ_Rogers(Orders,CoeffArr(:,i));
        
        if dX ==-100 && dY == -100 && dZ == -100
            exitFlagArr(i) = -1;
            endFlag = 1;
        else
            confinedShift = min(max([dX dY dZ],-1*[Q Q Q]),[Q Q Q]);
            PosArr(i,:) = PosArr(i,:) + confinedShift;
            if max(abs(PosArr(i,:))) > cropHalfSize
                exitFlagArr(i) = -2;
                endFlag = 1;
            elseif max(abs(confinedShift)) < Q
                if consecAccum == CritIter-1  %|| max(abs(confinedShift)) < 1e-5
                    TotPosArr(i,:) = PosArr(i,:) + maxXYZ(i,:); 
                    endFlag = 1;
                else
                    consecAccum = consecAccum + 1;
                end
            else
                consecAccum = 0;
            end
        end
    end
    fprintf(1,'peak %d, flag %d, consecAccum %d\n',i,exitFlagArr(i),consecAccum);
    toc
end
  
%% apply cutoff according to intensity
for i = 1:size(maxXYZ,1)
    if exitFlagArr(i) == 0
        goodAtomTotPos = TotPosArr(1:i-1,:);
        goodAtomTotPos = goodAtomTotPos(exitFlagArr(1:i-1)==0,:);
        Dist = pdist2(PosArr(i,:)+maxXYZ(i,:),goodAtomTotPos);
        if min(Dist) < minDist
            exitFlagArr(i) = -3;
        end
    end
end

%% classify the atom from non-atom
N                = 0;
exitFlagArr0     = exitFlagArr(:,1:end-N); 
TotPosArr0       = TotPosArr(1:end-N,:);
maxXYZ0          = maxXYZ(1:end-N,:);
atom_pos1        = TotPosArr0(exitFlagArr0==0,:); %
id               = find(exitFlagArr0==-3);
exitFlagArr0(id) = 0;
id0              = find(exitFlagArr0==-3);
minDist          = 0.9/Res;
for i = 1:size(maxXYZ0,1)
    if exitFlagArr0(i) == 0
        goodAtomTotPos = TotPosArr0(1:i-1,:);
        goodAtomTotPos = goodAtomTotPos(exitFlagArr0(1:i-1)==0,:);
        Dist = pdist2(PosArr(i,:)+maxXYZ0(i,:),goodAtomTotPos);
        if min(Dist) < minDist
            exitFlagArr0(i) = -3;
        end
    end
end
atom_pos0 = TotPosArr0(exitFlagArr0==0,:); 
atom_pos_tracing=atom_pos0;

%% initial classification
classify_info1 = struct('Num_species', 2,  'halfSize',  3,  'plothalfSize',  3, ...
      'O_Ratio', 1, 'SPHyn',  1,  'PLOT_YN',  1,  'separate_part',  40, 'lnorm',2);

[temp_model_ref, temp_atomtype_ref] = initial_class_kmean_sub_plot_ColorMap(FinalVol*10, atom_pos_tracing', classify_info1,3500);  
%% remove non-atoms
ID_type_nonAtoms                  = temp_atomtype_ref==1;
temp_model_C                      = atom_pos_tracing(:,:);
temp_model_C(ID_type_nonAtoms,:)  = [];
temp_atomtype_C                   = temp_atomtype_ref;
temp_atomtype_C(ID_type_nonAtoms) = [];
atom                              = temp_model_C'.*Res;
label                             = temp_atomtype_C;

save([pwd,'/output_data/atom_tracing_model.mat'],'atom','label')














