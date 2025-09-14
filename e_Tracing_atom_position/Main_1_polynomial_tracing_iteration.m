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
FinalVol = importdata([pwd,'/input_data/tBLG_reconstruction_volume.mat']);
FinalVol_org = FinalVol;
% padded the reconstruction with zero 
FinalVol_org = My_paddzero(FinalVol_org,size(FinalVol_org)+20,'double');


for iter = 1:21
    
    Dsetvol     = gather(FinalVol);                                        % the Dsetvol is 1646*1442*80 3D reconstruction object
    
    % Th: intensity threshold for the local maxima pixel
    % local maxima with intensity less than this value will not be traced
    % because they are way too weak to become actual atoms
    % MaxIter: the maximum iteration number for each atom
    % After MaxIter times iteration the engine will get out from the loop
    % CritIter: the maximum iteration before to check whether the local maxima is a good position for an atom
    % Res: the actual distance for each pixel
    % minDist: the minimum distance constraint in the unit of pixel
    % SearchiRad: the search radius for localizing an atom
    MaxIter = 14; CritIter = 7; interpol = 1;                     Th = 1.0e-2; 
    Res     = 0.19/2/interpol;   minDist = ceil(1.35/Res); SearchRad = 3;   
    
    % upsampling the reconstruction matrix by 3*3*3 by linear interpolation
    % better to run at super conputer since the size of the interpolated
    % volume will be larger than 16G
    xx = (1:size(Dsetvol,1)) - round((size(Dsetvol,1)+1)/2);               % xx yy zz is the vector around 0 the size of xx yy zz is (1, 3*N)
    yy = (1:size(Dsetvol,2)) - round((size(Dsetvol,2)+1)/2);
    zz = (1:size(Dsetvol,3)) - round((size(Dsetvol,3)+1)/2);
    
    xxi = ((interpol*xx(1)):(xx(end)*interpol))/interpol;                  % the value range of xxi yyi zzi are equal to xx yy zz
    yyi = ((interpol*yy(1)):(yy(end)*interpol))/interpol;                  % the size of xxi yyi zzi is (1, 3*N)
    zzi = ((interpol*zz(1)):(zz(end)*interpol))/interpol;
    
    % the size of  xxi yyi zzi is (1, 3*N-3), represent the sub pixel
    xxi = xxi(interpol:end); yyi = yyi(interpol:end); zzi = zzi(interpol:end); 
    [Y,X,Z]     = meshgrid(yy,xx,zz);                                      % generate the grid for yy xx zz(pixel) and yyi xxi zzi(sub pixel)
    [Yi,Xi,Zi]  = meshgrid(yyi,xxi,zzi);
    clear yyi  xxi zzi
    
    % upsampling the reconstruction matrix by 3*3*3 by linear interpolation
    Dsetvol     = interp3(Y,X,Z,Dsetvol,Yi,Xi,Zi,'spline',0);
    
    % padded the reconstruction with zero 
    FinalVol    = My_paddzero(Dsetvol,size(Dsetvol)+20,'double');
    
    % get polynomial power array
    fitCoeff = [];
    for tempi=0:4
        for j=0:4
            for k=0:4
                if tempi+j+k <= 4                
                    if max([tempi j k]) == 4
                        fitCoeff(end+1,:) = [tempi j k -1]; 
                    else                
                        fitCoeff(end+1,:) = [tempi j k 0]; 
                    end
                end
            end
        end
    end
    
    % get the local maxima from the reconstruction volume
    se = strel3d(3);
    dilatedBW   = imdilate(FinalVol,se);                                   % compare each pixel with it neighborhood
    maxPos      = find(FinalVol==dilatedBW & FinalVol>Th);                 % find the index for the local maxima pixel
    maxVals     = FinalVol(maxPos);                                        % obtain the value of each local maxima pixel 
    [~,sortInd] = sort(maxVals,'descend');                                 % sort the maxVals from large to small the first element is large
    maxNum      = min(100000,length(sortInd));                             % obtain the maximum number of potential position of atom
    maxPos      = maxPos(sortInd(1:maxNum));                               % update the local maxima pixel index with the maximum number
    
    % maxXYZ contain the (x,y,z) coordination position of the local maxima
    maxXYZ = zeros(length(maxPos),3);
    for tempi=1:length(maxPos)
        [xx,yy,zz] = ind2sub(size(FinalVol),maxPos(tempi));    
        if zz > 11 && zz < (size(FinalVol,3)-11)
            maxXYZ(tempi,:) = [xx yy zz];
        end
    end
    
    maxXYZ( ~any(maxXYZ,2), : ) = [];
    fprintf(1,'numpeak = %d \n',size(maxXYZ,1));
    clear Xi Yi Zi Dsetvol dilatedBW 
  
    % initialize the parameters
    Q = 0.5; Alpha = 5;                                                    % Q is the maximum sub pixel shift for each iteration in the loop % Alpha is the used to adjust the variance of the GaussWeight
    cropHalfSize   = SearchRad;                                            % the small region for polynomial Rogers fitting 7*7*7
    % obtain the coordination grid for the small region
    [X,Y,Z]      = ndgrid(-cropHalfSize:cropHalfSize,-cropHalfSize:cropHalfSize,-cropHalfSize:cropHalfSize);
    % obtain the index for the polynomial Roger fitting (now the region is a circle for out dataset we might have to change it into ellipse)
    SphereInd = find(X.^2+Y.^2+Z.^2 <=(SearchRad)^2);                      % obtiain the coordination grid for the polynomial Roger fitting
    % obtiain the coordination grid for the polynomial Roger fitting    
    XYZdata.X = X(SphereInd); XYZdata.Y = Y(SphereInd); XYZdata.Z = Z(SphereInd);               
    Orders    = fitCoeff(:,1:3);                                           % get polynomial power order
    PosArr    = zeros(size(maxXYZ));
    TotPosArr = zeros(size(maxXYZ));
    goodAtomTotPos = TotPosArr(1:0,:);
    % each Flag represent a state for an atom
    % -0: the atom position can be a correct position
    % -1: the fitting is wrong can only find the minimun not maximum
    % -2: the fitting maximum is beyond the searching small region
    % -3: the atom position is too close to the existed position 
    % -4: the iteration time is beyond the maximum iteration
    exitFlagArr  = zeros(1, size(maxXYZ,1));
    CoeffArr     = repmat(fitCoeff(:,4),[1 size(maxXYZ,1)]);
    % perform the main tracing loop
    countgoodatom = 0;

    for tempi =  1 :  +size(maxXYZ,1)
    
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
              exitFlagArr(tempi) = -4;
              endFlag = 1;
            end
            
            % obtain the grid position of the 7*7*7 small region for position searching
            cropXind = maxXYZ(tempi,1) + (-cropHalfSize:cropHalfSize);
            cropYind = maxXYZ(tempi,2) + (-cropHalfSize:cropHalfSize);
            cropZind = maxXYZ(tempi,3) + (-cropHalfSize:cropHalfSize);
            
            % obtain the value of the small region
            cropVol = FinalVol(cropXind,cropYind,cropZind);
            
            % obtian the fitting center position
            Pos = PosArr(tempi,:);
            % generate the GaussWeight for both expermential data and Roger fiting polynomial function
            GaussWeight = exp(-1*Alpha*( (X(SphereInd)-Pos(1)).^2 + (Y(SphereInd)-Pos(2)).^2 + (Z(SphereInd)-Pos(3)).^2) / (cropHalfSize)^2 );
            
            % define the fitting function
            fun = @(p,xdata) calculate_3D_polynomial_Rogers(xdata.X,xdata.Y,xdata.Z,Pos,Orders,p).*GaussWeight;
    
            opts = optimset('Display','off');
            
            % the last square fitting
            [p1,fminres1] = lsqcurvefit(fun,CoeffArr(:,tempi),XYZdata,cropVol(SphereInd).*GaussWeight,[],[],opts);
            % obtain the coefficient for current iteration
            CoeffArr(:,tempi) = p1;
            % calculate the local maximum
            [dX, dY, dZ] = calc_dX_dY_dZ_Rogers(Orders,CoeffArr(:,tempi));
    
            if dX ==-100 && dY == -100 && dZ == -100                           % only have local minimum which is wrong
                exitFlagArr(tempi) = -1;
                endFlag = 1;
            else
                % calculate the sub pixel shift value
                maxedShift = max([dX dY dZ],-1*[Q Q Q]);
                minedShift = min(maxedShift,[Q Q Q]);
                % shift the center
                PosArr(tempi,:) = PosArr(tempi,:) + minedShift;
                % center position beyond the small region
                if max(abs(PosArr(tempi,:))) > cropHalfSize
                    exitFlagArr(tempi) = -2;
                    endFlag = 1;
                % sub pixel shift is smaller than Q
                elseif max(abs(minedShift)) < Q
                    if consecAccum == CritIter-1                           % begin to judge the position is good or not
                        % the distance between current position with old good position
                        Dist = sqrt(sum((goodAtomTotPos - repmat(PosArr(tempi,:)+maxXYZ(tempi,:),[size(goodAtomTotPos,1) 1])).^2,2));
                        % the distance is too small 
                        if min(Dist) < minDist
                          exitFlagArr(tempi) = -3;
                        else
                          countgoodatom = countgoodatom+1;
                          TotPosArr(tempi,:) = PosArr(tempi,:) + maxXYZ(tempi,:);                       % store the good atom position
                          goodAtomTotPos(countgoodatom,:) = PosArr(tempi,:) + maxXYZ(tempi,:);          % store the good atom position
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
        fprintf(1,'peak %d, flag %d \n',tempi,exitFlagArr(tempi));
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
        %% save the final result
        FinalVol  = FinalVol_org; 
        atom(1,:) = atompos(1,:)-10;
        atom(2,:) = atompos(2,:)-10;
        atom(3,:) = atompos(3,:)-10;
        save([pwd,'/output_data/atom_tracing_model_total.mat'],'FinalVol','atom');
    else  
        atompos = goodAtomTotPos(:,:)';
    end
    vol_size = round(1/Res);TraceVol = Sub_Traced_Atom(volume,atompos,vol_size,iter);
    FinalVol = TraceVol(11:size(TraceVol,1)-10,11:size(TraceVol,2)-10,11:size(TraceVol,3)-10);                                                  
end

%% load the manual tracing result
clear;clc;  
% read in files: reconstruction volume
FinalVol = importdata([pwd,'/input_data/tBLG_reconstruction_volume.mat']);
FinalVol_org = FinalVol;
% padded the reconstruction with zero 
FinalVol_org = My_paddzero(FinalVol_org,size(FinalVol_org)+20,'double');
Res          = 0.19/2;
atom_auto    = importdata('atom_tracing_model_total.mat');atom_auto = atom_auto.atom; 
atom_manual  = importdata('atom_tracing_manual_model.mat');atom_manual = atom_manual.atom;
[l,w,h]=size(FinalVol_org(11:size(FinalVol_org,1)-10,11:size(FinalVol_org,2)-10,11:size(FinalVol_org,3)-10));
CentralPixels_Origin=([l,w,h]+1)/2;
atom_manual  = atom_manual./Res;
atom_manual  = atom_manual + CentralPixels_Origin';
% label out the manual tracing atom
manualLabel = ones(1,size(atom_manual,2));
for tempi=1:size(atom_manual,2)
    dif=(atom_auto-atom_manual(:,tempi));                                  % calculate all difference from the first result with tempi'th position in the second result (angstrom)
    dis=sqrt(sum(dif.^2,1));                                               % calculate the distance (angstrom)
    [dis,ind]=min(dis);                                                    % obatin the minimum distance and corresponding index
    if dis <= 1e-2                                                         % if the minimum distance is smaller than the threshold then store the information (here we set the minimum distance around 0.1 pm)
        manualLabel(tempi) = 1;
    else
        manualLabel(tempi) = 2;
    end
end
save([pwd,'/output_data/atom_tracing_manual_label.mat'],'manualLabel');



