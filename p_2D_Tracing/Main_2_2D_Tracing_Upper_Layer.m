%% Main Polynomial Tracing
% code to perform polynomial tracing on reconstruction volume
% reconstruction volume should be upsampled by 3*3*3 linear interpolation
% each local maximum is fitted with 9*9*9 voxels (3*3*3 before interpolation)
% by 4th order polynomial equation to get the position

clear;clc;
addpath([pwd,'/src/']);
addpath([pwd,'/input_data/']);
addpath([pwd,'/output_data/']);
% add the path to load reconstruction volume, you can comment it and move
% read in files: reconstruction volume
FinalVol = importdata([pwd,'/output_data/tBLG_reconstruction_volume.mat']);
Slice_org = FinalVol(:,:,21:21);
% padded the reconstruction with zero 
Slice_org = My_paddzero(Slice_org,size(Slice_org)+20,'double');
Slice     = Slice_org;


for iter = 1:21
    Res          = 0.19/2;
    minDist      = ceil(1.35/Res);
    se           = strel('square',3);
    dilatedBW    = imdilate(Slice, se);
    maxPos       = find( Slice == dilatedBW & Slice > 0.005);
    
    maxXYZ = zeros(length(maxPos),2);
    for i=1:length(maxPos)
        [xx,yy] = ind2sub(size(Slice),maxPos(i));
        maxXYZ(i,:) = [xx yy];
    end
    clear dilatedBW
    
    % Trace the atom from projection by 2D Gaussian
    cropHalfSize = 3;
    [bx,by]=ndgrid(-cropHalfSize:cropHalfSize,-cropHalfSize:cropHalfSize);
    BoxCoordinates.y=by;
    BoxCoordinates.x=bx;
    tracedPosition = zeros(size(maxXYZ));
    parfor i=1:size(maxXYZ,1)
        tic
        cropXind = maxXYZ(i,1) + (-cropHalfSize:cropHalfSize);
        cropYind = maxXYZ(i,2) + (-cropHalfSize:cropHalfSize);
        cropProj = Slice(cropXind,cropYind);
        fit_param_init  = [min(cropProj(:)), max(cropProj(:))-min(cropProj(:)), 0, 0, cropHalfSize/3, cropHalfSize/3, 0];
        fixed           = [0 ,  0,  0,  0,  0,    0,   0];
        lb              = [0 ,  0, -3, -3,  0,    0, -pi];
        ub              = [10, 10,  3,  3, inf, inf,  pi];
        [fit_resultI,~,~] = fit_gauss2D_PD(fit_param_init, BoxCoordinates, cropProj, fixed, lb, ub);
        tracedPosition(i,:) =  fit_resultI([4 3]);
        toc
    end
    tracedPosition = tracedPosition+maxXYZ;
    tracedPosition = tracedPosition';
    
    intensity = zeros(1,size(tracedPosition,2));
    parfor i=1:size(intensity,2)
        y = round(tracedPosition(1,i));
        x = round(tracedPosition(2,i));
        intensity(i) = sum(sum(Slice(y-1:y+1,x-1:x+1)));
    end
    [val,sortInd] = sort(intensity,'descend');
    tracedPosition = tracedPosition(:,sortInd);

    intensity = zeros(1,size(tracedPosition,2));
    parfor i=1:size(intensity,2)
        y = round(tracedPosition(1,i));
        x = round(tracedPosition(2,i));
        intensity(i) = Slice_org(y,x);
    end
    tracedPosition = tracedPosition(:,intensity>0.007);

    
    index = zeros(1,size(tracedPosition,2));
    for i=1:size(index,2)
        if index(i) == 0
            goodAtomTotPos = tracedPosition(:,1:i-1);
            goodAtomTotPos = goodAtomTotPos(:,index(1:i-1)==0);
            Dist = pdist2(tracedPosition(:,i)',goodAtomTotPos');
            if min(Dist)<minDist
                index(i) = 1;
            end
        end
    end
    tracedPosition = tracedPosition(:,index==0);
    
    if iter == 19
        atompos_1 = tracedPosition;
        atompos = atompos_1';
    elseif iter == 20
        atompos_2 = tracedPosition;
        atompos = [atompos_1,atompos_2];
    elseif iter == 21
        atompos_3 = tracedPosition;
        atompos = [atompos_1,atompos_2,atompos_3];     
        %% save the final result
        Slice  = Slice_org; 
        atom(1,:) = atompos(1,:)-10;
        atom(2,:) = atompos(2,:)-10;
        save([pwd,'/output_data/atom_tracing_model_Upper_Layer.mat'],'Slice','atom');
    else  
        atompos        = tracedPosition;
    end
    
    vol_size = round(1/Res);TraceSlice = Sub_Traced_Atom_2D(Slice_org,atompos,vol_size,1);
    Slice = TraceSlice;
end

load('atom_tracing_model_Upper_Layer.mat')
load('manual_atom_Upper_Layer.mat')
atom_manual = calculate_model_difference(atom,atom_manual,0.3/Res);
atom_manual = round(atom_manual);
index = zeros(1,size(atom_manual,2));
for tempi = 1:size(atom_manual,2)
    dif=(atom-atom_manual(:,tempi));
    dis=sqrt(sum(dif.^2,1));
    [val,ind] = min(dis);
    if val < 0.4/Res
        atom_manual(:,tempi) = atom(:,ind);
        index(tempi) = 1;
    end
end

tracedPosition = zeros(size(atom_manual));
parfor i=1:size(atom_manual,1)
    if index == 0
        tic
        cropXind = atom_manual(1,i) + (-cropHalfSize:cropHalfSize);
        cropYind = atom_manual(2,i) + (-cropHalfSize:cropHalfSize);
        cropProj = Slice(cropXind,cropYind);
        fit_param_init  = [min(cropProj(:)), max(cropProj(:))-min(cropProj(:)), 0, 0, cropHalfSize/3, cropHalfSize/3, 0];
        fixed           = [0 ,  0,  0,  0,  0,    0,   0];
        lb              = [0 ,  0, -3, -3,  0,    0, -pi];
        ub              = [10, 10,  3,  3, inf, inf,  pi];
        [fit_resultI,~,~] = fit_gauss2D_PD(fit_param_init, BoxCoordinates, cropProj, fixed, lb, ub);
        tracedPosition(:,i) =  fit_resultI([4 3])';
        toc
    end
end
atom = atom_manual + tracedPosition;
save([pwd,'/output_data/atom_tracing_model_Upper_Layer.mat'],'Slice','atom');

function atompos_2 = calculate_model_difference(atompos_1,atompos_2,atomDist)
Mag_shift=1;
shift=[0 0]';                                                            % the shift used to do the subpixel shift to align the atom position
while Mag_shift>1e-4                                                       % if the total difference is larger than 1e-6 (angstrom) then continue to do the alignment
    atompos_2 =  atompos_2 - shift;                                        % substruct the difference (i.e. do the shift to align the atom)
    difarr = [];                                                           % array for store the difference between two tracing result
    disarr = [];                                                           % array for store the distance between two tracing result
    count_arr1 = [];                                                       % store the index for the first atom result, which is match with the second atom tracing result
    count_arr2 = [];                                                       % store the index for the second atom result, which is match with the first atom tracing result
    for i=1:size(atompos_1,2)
        dif=(atompos_2-atompos_1(:,i));                                    % calculate all difference from the first result with i'th position in the second result (angstrom)
        dis=sqrt(sum(dif.^2,1));                                           % calculate the distance (angstrom)
        [dis,ind]=min(dis);                                                % obatin the minimum distance and corresponding index
        if dis <= atomDist                                                 % if the minimum distance is smaller than the threshold then store the information
            difarr=[difarr dif(:,ind)];
            disarr = [disarr dis];
            count_arr1 = [count_arr1,ind];
            count_arr2 = [count_arr2,i];
        end
    end
    shift=mean(difarr,2);                                                  % calculate the mean values of the difference in x-y-z axis, which is also the shift value for the next iteration
    Mag_shift=sum(abs(mean(difarr,2)));                                    % calculate the total difference
end
end



