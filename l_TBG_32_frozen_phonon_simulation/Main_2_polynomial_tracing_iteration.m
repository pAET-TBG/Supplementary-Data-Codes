%% Main Polynomial Tracing
% code to perform polynomial tracing on reconstruction volume
% reconstruction volume should be upsampled by 3*3*3 linear interpolation
% each local maximum is fitted with 9*9*9 voxels (3*3*3 before interpolation)
% by 4th order polynomial equation to get the position

clear;clc;
addpath('src/');
% add the path to load reconstruction volume, you can comment it and move
% the reconstruction into input folder
FinalVol = importdata([pwd,'/output_data/tBLG_reconstruction_volume.mat']);
FinalVol = FinalVol(20:end-20,20:end-20,:);
FinalVol_org = FinalVol;
% padded the reconstruction with zero 
FinalVol_org = My_paddzero(FinalVol_org,size(FinalVol_org)+20,'double');
% read in files: reconstruction volume

for iter = 1:23
    
    Dsetvol     = gather(FinalVol);                                                    % the Dsetvol is 300*300*300 3D reconstruction object
    
    % Th: intensity threshold for the local maxima pixel
    % local maxima with intensity less than this value will not be traced
    % because they are way too weak to become actual atoms
    % MaxIter: the maximum iteration number for each atom
    % After MaxIter times iteration the engine will get out from the loop
    % CritIter: the maximum iteration before to check whether the local maxima is a good position for an atom
    % Res: the actual distance for each pixel
    % minDist: the minimum distance constraint in the unit of pixel
    % SearchiRad: the search radius for localizing an atom
    MaxIter = 14;       CritIter = 7;       interpol = 1;  Th = 0.5e-2; 
    Res     = 0.19/2/interpol;   minDist  = round(1.4/Res);     
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
    clear yyi  xxi zzi xx yy zz
    
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
        fprintf(1,'peak %d, flag %d \n',i,exitFlagArr(i));
        fprintf(['number of good atom: ', num2str(countgoodatom),'\n']);
    end

    %% z-direction constraint
    Index = (goodAtomTotPos(:,3) >= 63 & goodAtomTotPos(:,3) <= 73) | (goodAtomTotPos(:,3) >= 27 & goodAtomTotPos(:,3) <= 37);
    goodAtomTotPos = goodAtomTotPos(Index,:);

    %% substract the traced atom from volume and prepare FinalVol for the next iteration
    volume = FinalVol_org;      
    if iter == 19
        atompos_1 = goodAtomTotPos;
        atompos = atompos_1';
    elseif iter == 20
        atompos_2 = goodAtomTotPos;
        atompos = [atompos_1',atompos_2'];
    elseif iter == 21
        atompos_3 = goodAtomTotPos;
        atompos = [atompos_1',atompos_2',atompos_3']; 
    elseif iter == 22
        atompos_4 = goodAtomTotPos;
        atompos = [atompos_1',atompos_2',atompos_3',atompos_4'];
    elseif iter == 23
        atompos_5 = goodAtomTotPos;
        atompos = [atompos_1',atompos_2',atompos_3',atompos_4',atompos_5'];
        %% save the final result
        FinalVol  = FinalVol_org; 
        atom(1,:) = atompos(1,:)-10;
        atom(2,:) = atompos(2,:)-10;
        atom(3,:) = atompos(3,:)-10;
        save([pwd,'/output_data/atom_tracing_model_total.mat'],'FinalVol','atom');
    else  
        atompos = goodAtomTotPos(:,:)';
    end

    %% pick up the average volume for atom
    if iter <= 6
        % with small iteration we want to find all possible position 
        sigma_filter = 5;
    elseif iter >= 16
        % with large iteration we want the tracing result converge
        sigma_filter = 6;
    else
        % gradually enlarge the atom size to elimilate the wrong position
        sigma_filter = 0.1*iter + 4.4;
    end
    vol_size = round(1/Res);
    atom_vol = zeros(2*vol_size+1,2*vol_size+1,2*vol_size+1);
    filter = create3DGeneralizedGaussianFilter(vol_size,sigma_filter,2);
    filter = filter./max(filter(:)).*1;
    count = 0;
    for tempi = 1:size(atompos,2)
        if ((round(atompos(2,tempi))-vol_size)>0) && ((round(atompos(1,tempi))-vol_size)>0) && ((round(atompos(3,tempi))-vol_size)>0) ...
           && (round(atompos(2,tempi))+vol_size<size(volume,2)) && (round(atompos(1,tempi))+vol_size<size(volume,1)) && (round(atompos(3,tempi))+vol_size<size(volume,3))
            x = round(atompos(2,tempi))-vol_size:round(atompos(2,tempi))+vol_size;
            y = round(atompos(1,tempi))-vol_size:round(atompos(1,tempi))+vol_size;
            z = round(atompos(3,tempi))-vol_size:round(atompos(3,tempi))+vol_size;
            atom_vol = atom_vol + volume(y,x,z).*filter;
            count = count + 1;
        end
    end
    atom_vol = atom_vol./count;

    %% substact the traced atom
    TraceVol = volume;
    for tempi = 1:size(atompos,2)
        if ((round(atompos(2,tempi))-vol_size)>0) && ((round(atompos(1,tempi))-vol_size)>0) && ((round(atompos(3,tempi))-vol_size)>0) ...
           && (round(atompos(2,tempi))+vol_size<size(volume,2)) && (round(atompos(1,tempi))+vol_size<size(volume,1)) && (round(atompos(3,tempi))+vol_size<size(volume,3))
            x = round(atompos(2,tempi))-vol_size:round(atompos(2,tempi))+vol_size;
            y = round(atompos(1,tempi))-vol_size:round(atompos(1,tempi))+vol_size;
            z = round(atompos(3,tempi))-vol_size:round(atompos(3,tempi))+vol_size;
            TraceVol(y,x,z) = TraceVol(y,x,z) - atom_vol;
        end
    end

    %% convolution to generate new volume
    % generate a 3D Gaussian filter to smooth the volume %% before 6 I use 5
    filter = atom_vol.*create3DGeneralizedGaussianFilter(vol_size,sigma_filter,2);filter = filter./sum(filter(:));
    TraceVol = convn(TraceVol,filter,'same');
    TraceVol = TraceVol(11:size(TraceVol,1)-10,11:size(TraceVol,2)-10,11:size(TraceVol,3)-10);                                                  
    FinalVol = TraceVol;

end

%% matched atom result
%% load the manual tracing result
clear;clc;
addpath([pwd,'/src/']);
addpath([pwd,'/input_data/']);
addpath([pwd,'/output_data/']);
% read in files: reconstruction volume
FinalVol = importdata([pwd,'/output_data/tBLG_reconstruction_volume.mat']);
FinalVol = FinalVol(20:end-20,20:end-20,:);
FinalVol_org = FinalVol;
% padded the reconstruction with zero 
FinalVol_org = My_paddzero(FinalVol_org,size(FinalVol_org)+20,'double');
Res          = 0.19/2;
atom_auto    = importdata('atom_tracing_model_total.mat');
atom_auto    = atom_auto.atom; 
atom_manual  = importdata('atom_tracing_manual_model_total.mat');
label        = atom_manual.label;
atom_manual  = atom_manual.atom;
[l,w,h]=size(FinalVol_org(11:size(FinalVol_org,1)-10,11:size(FinalVol_org,2)-10,11:size(FinalVol_org,3)-10));
CentralPixels_Origin=([l,w,h]+1)/2;
atom_auto    = atom_auto - CentralPixels_Origin';
atom_auto    = atom_auto.*Res;

%% match the atom
atomDist = 0.01;
atom = atom_manual;
count = 0;
for i=1:size(atom_manual,2)
    dif=(atom_auto-atom_manual(:,i));          
    dis=sqrt(sum(dif.^2,1));   
    [dis,ind]=min(dis);    
    if dis <= atomDist    
        count = count+1;
        atom(:,i) = atom_auto(:,ind);
    end
end
save([pwd,'/output_data/atom_tracing_manual_model.mat'],'atom','label')

function filter = create3DGeneralizedGaussianFilter(size, sigma, beta)
    % size is the length of each dimension of the filter
    % sigma is the standard deviation of the Gaussian distribution

    % Create a grid of coordinates
    [x, y, z] = ndgrid(-size:size);

    % Calculate the Gaussian function
    filter = exp(-abs((x.^2 + y.^2 + z.^2) / (beta * sigma^2)).^beta);

    % Normalize the filter so that the sum of all elements is 1
    filter = filter ./ sum(filter(:));
end
