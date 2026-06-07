clear;
clc;

%% load data
addpath([pwd,'\src\'])
addpath([pwd,'\input_data\'])
addpath([pwd,'\output_data\'])

projs     = importdata('Projections_Initial.mat');
angles    = importdata('Angles.mat');
Vol       = importdata('Volume_Atom_Tracing.mat');

%% prepare the volume
Vol = My_paddzero(Vol,[2000,2000,104]);
Vol = imtranslate(Vol,[-18,10,0]);

%% rotation setting
rotation       = 'ZYX';  % Euler angles setting ZYZ
dtype          = 'single';
angles_refined      = cast(angles,dtype);
Num_pj = size(angles_refined,1);

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

% extract size of projections & num of projections
vec1 = matR(:,1);
vec2 = matR(:,2);
vec3 = matR(:,3);

%% rotation matrix
Rs = zeros(3,3,Num_pj, dtype);
for k = 1:Num_pj
    phi   = angles_refined(k,1);
    theta = angles_refined(k,2);
    psi   = angles_refined(k,3);
    
    % compute rotation matrix R w.r.t euler angles {phi,theta,psi}
    rotmat1 = MatrixQuaternionRot(vec1,phi);
    rotmat2 = MatrixQuaternionRot(vec2,theta);
    rotmat3 = MatrixQuaternionRot(vec3,psi);
    R =  single(rotmat1*rotmat2*rotmat3)';
    Rs(:,:,k) = R;
end

%% compute calculated projections
ref_projs = calculate3Dprojection_multiGPU(single(Vol), Rs);
ref_projs = ref_projs(750:1250, 510:1250,:);
projs_sm  = projs(750:1250, 510:1250,:);


%% normalize each projection
for tempi = 1:size(projs,3)
    tempProj   = projs_sm(:,:,tempi);
    
    projs(:,:,tempi)     = projs(:,:,tempi) - min(tempProj(:));
    ref_projs(:,:,tempi) = ref_projs(:,:,tempi) - min(min(ref_projs(:,:,tempi)));

    projs(:,:,tempi)     = projs(:,:,tempi) ./max(tempProj(:)-min(tempProj(:)));
    ref_projs(:,:,tempi) = ref_projs(:,:,tempi) ./max(max(ref_projs(:,:,tempi)));

    clear tempProj tempProj_d
end



%% cross-correlation alignment
for iter = 1:2
    projs_sm = projs(755:1245, 515:1245,:);
    for proj_index = 1:size(projs,3)
    
        %% cross correlation to match the experiment raw data
        [colS, rowS] = size(projs_sm(:,:,proj_index));
        [colB, rowB] = size(ref_projs(:,:,proj_index));
        corrMatrix = double(zeros(colB - colS + 1, rowB - rowS +1)); 
        for colPos = 1:colB - colS + 1 
            for rowPos = 1:rowB - rowS + 1        
                tempProj = ref_projs(colPos:colPos+colS-1, rowPos:rowPos+rowS-1,proj_index);
                corrMatrix(colPos, rowPos) = corr2(tempProj, projs_sm(:,:,proj_index));
            end
        end
        
        %% find the maximum position
        [maxNum1, tempR] = max(max(corrMatrix));
        [maxNum2, tempC] = max(corrMatrix(:,tempR));
        maxProj =  ref_projs(tempC:tempC+colS-1, tempR:tempR+rowS-1,proj_index);

        %% update the projection
        projs(:,:,proj_index) = imtranslate(projs(:,:,proj_index), ...
            [tempR-round(size(corrMatrix,2)+1)/2,tempC-round(size(corrMatrix,1)+1)/2],'cubic');        
    end
end
clear projs_sm  ref_projs

%% normalization the projection
% projs       = projs_org;
projs_sm    = projs(750:1250, 510:1250,:);
ref_projs   = calculate3Dprojection_multiGPU(single(Vol), Rs);
ref_projs   = ref_projs(750:1250, 510:1250,:);
ref_projs   = ref_projs./max(ref_projs(:));

%% match the intensity for each projection again
for tempi = 1:size(projs,3)
    % obtain the simulation and experiment projections
    tempProj_sim   = ref_projs(:,:,tempi);
    tempProj_exp   = projs_sm(:,:,tempi);

    % normalization
    tempProj_sim       = tempProj_sim       - min(tempProj_sim(:));
    tempProj_sim       = tempProj_sim       ./max(tempProj_sim(:));
    tempProj_exp       = tempProj_exp       - min(tempProj_exp(:));
    tempProj_exp       = tempProj_exp       ./max(tempProj_exp(:)); 
    
    % obtain the distribution of intensity
    [countsim, IntSim]     = imhist(tempProj_sim,1000);
    [countexp, IntExp]     = imhist(tempProj_exp,1000);

    % find the intensity that at list 50 pixels in the intensity range
    MatchNum     = max(countsim(:))./2;
    indicesSim   = find(countsim > MatchNum);
    indciesExp   = find(countexp > MatchNum);
    sinMin       = IntSim(indicesSim(1));
    sinMax       = IntSim(indicesSim(end));
    expMin       = IntSim(indciesExp(1));
    expMax       = IntSim(indciesExp(end));
    
    % match the intensity of experiment raw data with simulation data
    tempProj_exp     = (tempProj_exp - expMin)./expMax.*sinMax + sinMin;
    projs(:,:,tempi) = (projs(:,:,tempi) - expMin)./expMax.*sinMax + sinMin;
    projs(:,:,tempi) = projs(:,:,tempi) - mean(tempProj_exp(:)) + mean(tempProj_sim(:));
    
    % de-normalization for experiment raw data
    projs(:,:,tempi) = projs(:,:,tempi).*max(max(ref_projs(:,:,tempi) - min(min(ref_projs(:,:,tempi))))) ...
                       + min(min(ref_projs(:,:,tempi)));
    
    clear tempProj_exp  tempProj_sim
end

save([pwd,'\output_data\Projections.mat'],'projs');



