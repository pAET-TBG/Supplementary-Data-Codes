%% Main Polynomial Tracing
% code to perform polynomial tracing on reconstruction volume
% reconstruction volume should be upsampled by 3*3*3 linear interpolation
% each local maximum is fitted with 9*9*9 voxels (3*3*3 before interpolation)
% by 4th order polynomial equation to get the position

clear;clc;
addpath('src/');
% add the path to load reconstruction volume, you can comment it and move
% the reconstruction into input folder
FinalVol = importdata([pwd,'/output_data/Rec_Without_Probe_NVcenter_Multislice.mat']);
FinalVol = FinalVol(17:227,265:485,30:770);
FinalVol_1 = FinalVol(:,:,1:258);
FinalVol_2 = FinalVol(:,:,246:506);
FinalVol_3 = FinalVol(:,:,495:741);
clear FinalVol


for iter = 1:3
    if iter == 1
        Dsetvol = gather(FinalVol_1);                                      % the Dsetvol is 300*300*300 3D reconstruction object
    elseif iter == 2
        Dsetvol = gather(FinalVol_2);                                      % the Dsetvol is 300*300*300 3D reconstruction object
    elseif iter == 3
        Dsetvol = gather(FinalVol_3);                                      % the Dsetvol is 300*300*300 3D reconstruction object
    end
    % Th: intensity threshold for the local maxima pixel
    % local maxima with intensity less than this value will not be traced
    % because they are way too weak to become actual atoms
    % MaxIter: the maximum iteration number for each atom
    % After MaxIter times iteration the engine will get out from the loop
    % CritIter: the maximum iteration before to check whether the local maxima is a good position for an atom
    % Res: the actual distance for each pixel
    % minDist: the minimum distance constraint in the unit of pixel
    % SearchiRad: the search radius for localizing an atom
    MaxIter = 14;       CritIter = 7;       interpol = 3;  Th = 5.5e-3; 
    Res     = 0.2129/interpol;   minDist  = round(1/Res);     
    SearchRadxy = 3;
    SearchRadz  = 3;
    
    % upsampling the reconstruction matrix by 3*3*3 by linear interpolation
    % better to run at super conputer since the size of the interpolated
    % volume will be larger than 16G
    xx = (1:size(Dsetvol,1)) - round((size(Dsetvol,1)+1)/2);                   % xx yy zz is the vector around 0 the size of xx yy zz is (1, 3*N)
    yy = (1:size(Dsetvol,2)) - round((size(Dsetvol,2)+1)/2);
    zz = (1:size(Dsetvol,3)) - round((size(Dsetvol,3)+1)/2);
    
    xxi = ((interpol*xx(1)):(xx(end)*interpol))/interpol;                      % the value range of xxi yyi zzi are equal to xx yy zz
    yyi = ((interpol*yy(1)):(yy(end)*interpol))/interpol;                      % the size of xxi yyi zzi is (1, 3*N)
    zzi = ((interpol*zz(1)):(zz(end)*interpol))/interpol;
    
    xxi = xxi(interpol:end); yyi = yyi(interpol:end); zzi = zzi(interpol:end); % the size of  xxi yyi zzi is (1, 3*N-3), represent the sub pixel
    [Y,X,Z]     = meshgrid(yy,xx,zz);                                          % generate the grid for yy xx zz(pixel) and yyi xxi zzi(sub pixel)
    [Yi,Xi,Zi]  = meshgrid(yyi,xxi,zzi);
    clear yyi  xxi zzi
    
    % upsampling the reconstruction matrix by 3*3*3 by linear interpolation
    Dsetvol     = interp3(Y,X,Z,Dsetvol,Yi,Xi,Zi,'spline',0);
    
    % padded the reconstruction with zero 
    FinalVol    = My_paddzero(Dsetvol,size(Dsetvol)+20,'double');
    
    % get polynomial power array
    fitCoeff = [];
    for i=0:4
        for j=0:4
            for k=0:4
                if i+j+k <= 4                
                    if max([i j k]) == 4
                        fitCoeff(end+1,:) = [i j k -1];                %#ok<SAGROW>
                    else                
                        fitCoeff(end+1,:) = [i j k 0];                 %#ok<SAGROW>
                    end
                end
            end
        end
    end
    
    % get the local maxima from the reconstruction volume
    % se = strel3dEllpse(3,5);
    se = strel3d(3);
    dilatedBW   = imdilate(FinalVol,se);                                       % compare each pixel with it neighborhood
    maxPos      = find(FinalVol==dilatedBW & FinalVol>Th);                     % find the index for the local maxima pixel
    maxVals     = FinalVol(maxPos);                                            % obtain the value of each local maxima pixel 
    [~,sortInd] = sort(maxVals,'descend');                                     % sort the maxVals from large to small the first element is large
    maxNum      = min(100000,length(sortInd));                                 % obtain the maximum number of potential position of atom
    maxPos      = maxPos(sortInd(1:maxNum));                                   % update the local maxima pixel index with the maximum number
    
    % maxXYZ contain the (x,y,z) coordination position of the local maxima
    maxXYZ = zeros(length(maxPos),3);
    for i=1:length(maxPos)
        [xx,yy,zz] = ind2sub(size(FinalVol),maxPos(i));    
        if zz > 11 && zz < (size(FinalVol,3)-11)
            maxXYZ(i,:) = [xx yy zz];
        end
    end
    
    maxXYZ( ~any(maxXYZ,2), : ) = [];
    fprintf(1,'numpeak = %d \n',size(maxXYZ,1));
    clear Xi Yi Zi Dsetvol dilatedBW 
    
    % initialize the parameters
    Qxy = 0.5;  Qz = 0.5;   Alpha = 5;                                         % Q is the maximum sub pixel shift for each iteration in the loop % Alpha is the used to adjust the variance of the GaussWeight
    cropHalfSizexy = SearchRadxy;                                              % the small region for polynomial Rogers fitting 7*7*15
    cropHalfSizez  = SearchRadz;
    [X,Y,Z] = ndgrid(-cropHalfSizexy:cropHalfSizexy, ...                       % obtain the coordination grid for the small region
                     -cropHalfSizexy:cropHalfSizexy, ...
                     -cropHalfSizez:cropHalfSizez);
    SphereInd = find(X.^2+Y.^2+(Z./SearchRadz.*SearchRadxy).^2 <=(2.*SearchRadxy./3+SearchRadz./3)^2);           % obtain the index for the polynomial Roger fitting (now the region is a circle for out dataset we might have to change it into ellipse)
    XYZdata.X = X(SphereInd); XYZdata.Y = Y(SphereInd); XYZdata.Z = Z(SphereInd);                % obtiain the coordination grid for the polynomial Roger fitting
    Orders       = fitCoeff(:,1:3);                                                              % get polynomial power order
    PosArr = zeros(size(maxXYZ));
    TotPosArr = zeros(size(maxXYZ));
    goodAtomTotPos = TotPosArr(1:0,:);
    % each Flag represent a state for an atom
    % 0 : the atom position can be a correct position
    % -1: the fitting is wrong can only find the minimun not maximum
    % -2: the fitting maximum is beyond the searching small region
    % -3: the atom position is too close to the existed position 
    % -4: the iteration time is beyond the maximum iteration
    % exitFlagArr  = zeros(1, size(maxXYZ1,1)+size(maxXYZ2,1)+size(maxXYZ3,1)+size(maxXYZ4,1));
    % CoeffArr     = repmat(fitCoeff(:,4),[1 size(maxXYZ1,1)+size(maxXYZ2,1)+size(maxXYZ3,1)+size(maxXYZ4,1)]);
    exitFlagArr  = zeros(1, size(maxXYZ,1));
    CoeffArr     = repmat(fitCoeff(:,4),[1 size(maxXYZ,1)]);
    % perform the main tracing loop
    totalAtomNum = 0;
    countgoodatom = 0;
    
    tic
    for i = totalAtomNum+1 : totalAtomNum+size(maxXYZ,1)
        % if the endFlag is equal to 0 then we still in the loop to find the atom position    
        endFlag = 0;
        % count for the operate times when sub pixel shift is smaller than Q
        % when consecAccum == CritIter-1 begin to check whether the atom is good
        consecAccum = 0;
        % the iterNum should be smaller than MaxIter
        iterNum = 0;
        
        while ~endFlag    
            iterNum = iterNum + 1;
            if iterNum>MaxIter                                                 % if the iterNUm is larger than MaxIter then break from the loop
              exitFlagArr(i) = -4;
              endFlag = 1;
            end
            
            % obtain the grid position of the 7*7*7 small region for position searching
            cropXind = maxXYZ(i-totalAtomNum,1) + (-cropHalfSizexy:cropHalfSizexy);
            cropYind = maxXYZ(i-totalAtomNum,2) + (-cropHalfSizexy:cropHalfSizexy);
            cropZind = maxXYZ(i-totalAtomNum,3) + (-cropHalfSizez:cropHalfSizez);
            
            % obtain the value of the small region
            cropVol = FinalVol(cropXind,cropYind,cropZind);
            
            % obtian the fitting center position
            Pos = PosArr(i,:);
            % generate the GaussWeight for both expermential data and Roger fiting polynomial function
            GaussWeight = exp(-1*Alpha*( (X(SphereInd)-Pos(1)).^2 + (Y(SphereInd)-Pos(2)).^2 + (Z(SphereInd)-Pos(3)).^2) / (cropHalfSizexy)^2 );
            
            % define the fitting function
            fun = @(p,xdata) calculate_3D_polynomial_Rogers(xdata.X,xdata.Y,xdata.Z,Pos,Orders,p).*GaussWeight;
    
            opts = optimset('Display','off');
            
            % the last square fitting
            [p1,fminres1] = lsqcurvefit(fun,CoeffArr(:,i),XYZdata,cropVol(SphereInd).*GaussWeight,[],[],opts);
            % obtain the coefficient for current iteration
            CoeffArr(:,i) = p1;
            % calculate the local maximum
            [dX, dY, dZ] = calc_dX_dY_dZ_Rogers(Orders,CoeffArr(:,i));
    
            if dX ==-100 && dY == -100 && dZ == -100                           % only have local minimum which is wrong
                exitFlagArr(i) = -1;
                endFlag = 1;
            else
                % calculate the sub pixel shift value
                maxedShift = [max([dX dY],-1*[Qxy Qxy]),max(dZ,-1*Qz)];
                minedShift = [min(maxedShift(1:2),[Qxy Qxy]),min(maxedShift(3),Qz)];
                % shift the center
                PosArr(i,:) = PosArr(i,:) + minedShift;
                % center position beyond the small region
                if (max(abs(PosArr(i,1:2))) > cropHalfSizexy) || (abs(PosArr(i,3)) > cropHalfSizez)
                    exitFlagArr(i) = -2;
                    endFlag = 1;
                % sub pixel shift is smaller than Q
                elseif (max(abs(minedShift(1:2))) < Qxy) && (abs(minedShift(3)) < Qz)
                    if consecAccum == CritIter-1                               % begin to judge the position is good or not
                        if i == 1
                            goodAtomTotPos = TotPosArr(1:i-1,:);                   % obtain all good atom position
                            goodAtomTotPos = goodAtomTotPos(exitFlagArr(1:i-1)==0,:);
                        end
                        % the distance between current position with old good position
                        Dist = sqrt(sum(([goodAtomTotPos(:,1:2),goodAtomTotPos(:,3)] - ...
                                         [repmat(PosArr(i,1:2)+maxXYZ(i-totalAtomNum,1:2),[size(goodAtomTotPos,1) 1]), ...
                                          repmat((PosArr(i,3)+maxXYZ(i-totalAtomNum,3)),[size(goodAtomTotPos,1) 1])]).^2,2));
                        % the distance is too small 
                        if min(Dist) < minDist
                          exitFlagArr(i) = -3;
                        else
                          countgoodatom = countgoodatom+1;
                          TotPosArr(i,:) = PosArr(i,:) + maxXYZ(i-totalAtomNum,:);          % store the good atom position
                          goodAtomTotPos(countgoodatom,:) = PosArr(i,:) + maxXYZ(i-totalAtomNum,:);          % store the good atom position
                          cropFitVol = zeros(2*cropHalfSizexy+1,2*cropHalfSizexy+1,2*cropHalfSizez+1);
                          cropFitVol(SphereInd) = fun(CoeffArr(:,i),XYZdata);
                          clear cropFitVol
                        end
                        endFlag = 1;
                    else
                        consecAccum = consecAccum + 1;                         % count when sub pixel shift is smaller than Q
                    end
                else
                    consecAccum = 0;                                           % id the sub pixel shift is larger than Q then start from the beginning
                end
            end
        end
        if mod(i,100) == 0
            toc
        end
        fprintf(1,'peak %d, flag %d \n',i,exitFlagArr(i));
        fprintf(['number of good atom: ', num2str(countgoodatom),'\n']);
    end

    %% save this iteration traced atom0
    save([pwd,'/output_data/atom_tracing_NVC_Without_Probe_Multislice_',num2str(iter),'.mat'],'goodAtomTotPos','FinalVol','-v7.3');
    clear FinalVol
end

%% load the auto tracing result 
load('atom_tracing_NVC_Without_Probe_Multislice_1.mat');
atom1 = goodAtomTotPos';
load('atom_tracing_NVC_Without_Probe_Multislice_2.mat');
atom2 = goodAtomTotPos';
load('atom_tracing_NVC_Without_Probe_Multislice_3.mat');
atom3 = goodAtomTotPos';
atom2(3,:) = atom2(3,:) + (246-1)*interpol;
atom3(3,:) = atom3(3,:) + (495-1)*interpol;
% match the pixel size
atom1 = atom1./interpol;
atom2 = atom2./interpol;
atom3 = atom3./interpol;

%% finely align and match different region togather
atom1_align = atom1(:,atom1(3,:)<=254&atom1(3,:)>=250);
atom2_align = atom2(:,atom2(3,:)<=254&atom2(3,:)>=250);
atomDist  = 0.5;
Mag_shift = 1;
shift=[0 0 0]';                                                            % the shift used to do the subpixel shift to align the atom position
while Mag_shift> 1e-5                                                      % if the total difference is larger than 1e-6 (angstrom) then continue to do the alignment
    atom2 = atom2-shift;
    atom2_align = atom2_align  - shift;
    difarr = [];                                                           % array for store the difference between two tracing result
    for i=1:size(atom1_align,2)
        dif=(atom2_align-atom1_align(:,i));                                % calculate all difference from the first result with i'th position in the second result (angstrom)
        dis=sqrt(sum(dif.^2,1));                                           % calculate the distance (angstrom)
        [dis,ind]=min(dis);                                                % obatin the minimum distance and corresponding index
        if dis <= atomDist                                                 % if the minimum distance is smaller than the threshold then store the information
            difarr=[difarr dif(:,ind)];
        end
    end
    shift=mean(difarr,2);                                                  % calculate the mean values of the difference in x-y-z axis, which is also the shift value for the next iteration
    Mag_shift=sum(abs(mean(difarr,2)));                                    % calculate the total difference
end

atom2_align = atom2(:,atom2(3,:)<=503&atom2(3,:)>=498);
atom3_align = atom3(:,atom3(3,:)<=503&atom3(3,:)>=498);
atomDist  = 0.5;
Mag_shift = 1;
shift=[0 0 0]';                                                            % the shift used to do the subpixel shift to align the atom position
while Mag_shift> 1e-5                                                      % if the total difference is larger than 1e-6 (angstrom) then continue to do the alignment
    atom3 = atom3-shift;
    atom3_align = atom3_align  - shift;
    difarr = [];                                                           % array for store the difference between two tracing result                                              
    for i=1:size(atom1_align,2)
        dif=(atom3_align-atom2_align(:,i));                                % calculate all difference from the first result with i'th position in the second result (angstrom)
        dis=sqrt(sum(dif.^2,1));                                           % calculate the distance (angstrom)
        [dis,ind]=min(dis);                                                % obatin the minimum distance and corresponding index
        if dis <= atomDist                                                 % if the minimum distance is smaller than the threshold then store the information
            difarr=[difarr dif(:,ind)];
        end
    end
    shift=mean(difarr,2);                                                  % calculate the mean values of the difference in x-y-z axis, which is also the shift value for the next iteration
    Mag_shift=sum(abs(mean(difarr,2)));                                    % calculate the total difference
end

atom = [atom1 atom2 atom3];
save([pwd,'/output_data/atom_tracing_NVC_Without_Probe_Multislice_total.mat'],'atom','-v7.3');