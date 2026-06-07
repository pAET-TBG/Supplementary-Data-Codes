clear;
clc;

%% load data
addpath([pwd,'\src\'])
addpath([pwd,'\input_data\'])
addpath([pwd,'\output_data\'])

projs     = importdata('Projections_Initial.mat');
angles    = importdata('Angles_Initial.mat');
Vol       = importdata('Volume_Atom_Tracing.mat');

%% set the projections' out edge to be zero
projs(1:291,:,1) = 0;projs(:,1:307,1) = 0;projs(1903:2000,:,1) = 0;projs(:,1901:2000,1) = 0;
projs(1:206,:,2) = 0;projs(:,1:284,2) = 0;projs(1820:2000,:,2) = 0;projs(:,1895:2000,2) = 0;
projs(1:172,:,3) = 0;projs(:,1:372,3) = 0;projs(1783:2000,:,3) = 0;projs(:,1979:2000,3) = 0;
projs(1:165,:,4) = 0;projs(:,1:211,4) = 0;projs(1774:2000,:,4) = 0;projs(:,1822:2000,4) = 0;
projs(1:208,:,5) = 0;projs(:,1:311,5) = 0;projs(1817:2000,:,5) = 0;projs(:,1922:2000,5) = 0;
projs(1:131,:,6) = 0;projs(:,1:370,6) = 0;projs(1729:2000,:,6) = 0;projs(:,1961:2000,6) = 0;
projs(1:147,:,7) = 0;projs(:,1:281,7) = 0;projs(1758:2000,:,7) = 0;projs(:,1893:2000,7) = 0;
projs(1:243,:,8) = 0;projs(:,1:279,8) = 0;projs(1851:2000,:,8) = 0;projs(:,1886:2000,8) = 0;
projs(1:115,:,9) = 0;projs(:,1:200,9) = 0;projs(1731:2000,:,9) = 0;projs(:,1809:2000,9) = 0;
projs(1:358,:,10) = 0;projs(:,1:295,10) = 0;projs(1970:2000,:,10) = 0;projs(:,1911:2000,10) = 0;
projs(1:229,:,11) = 0;projs(:,1:143,11) = 0;projs(1838:2000,:,11) = 0;projs(:,1752:2000,11) = 0;
projs(1:151,:,12) = 0;projs(:,1:261,12) = 0;projs(1763:2000,:,12) = 0;projs(:,1868:2000,12) = 0;
projs(1:214,:,13) = 0;projs(:,1:285,13) = 0;projs(1824:2000,:,13) = 0;projs(:,1893:2000,13) = 0;
    

%% roughly align the stable marker and recenter the projection and reference volume
Vol = My_paddzero(Vol,[2000 2000 104]);
Vol = imtranslate(Vol,[-17,7,0]);

% [-2:1:2] -> [-1:0.5:1] -> [-0.5:0.1:0.5]-> [-0.1:0.05:0.1]
% [-0.2:0.1:0.2] -> [-0.1:0.1:0.1] -> [-0.1:0.1:0.1]-> [-0.1:0.05:0.1]
count = 1;
scanRange_xy{1} = -1:1:1;scanRange_xy{2} = -0.5:0.5:0.5;scanRange_xy{3} = -0.25:0.25:0.25;scanRange_xy{4} = -0.1:0.1:0.1;scanRange_xy{5} = -0.05:0.05:0.05; 
scanRange_z{1} = -0.1:0.1:0.1;scanRange_z{2} = -0.1:0.1:0.1;scanRange_z{3} = -0.1:0.1:0.1;scanRange_z{4} = -0.1:0.1:0.1;scanRange_z{5} = -0.05:0.05:0.05;   

while count ~= 6
    for index = 1:size(projs,3)
        
        angle_center = angles(index,:);
        xAxis_angles = scanRange_xy{count}; 
        yAxis_angles = scanRange_xy{count};
        zAxis_angles = scanRange_z{count};
        angle_Matrix = zeros(size(xAxis_angles,2),size(yAxis_angles,2),size(zAxis_angles,2));
        
        for tempx = 1:size(xAxis_angles,2)
            for tempy = 1:size(yAxis_angles,2)
                for tempz = 1:size(zAxis_angles,2)
                    tic
                    angle(1,1) = angle_center(1)+zAxis_angles(tempz);
                    angle(1,2) = angle_center(2)+xAxis_angles(tempx);
                    angle(1,3) = angle_center(3)+yAxis_angles(tempy);
        
                    %% rotation setting
                    rotation       = 'ZYX';  % Euler angles setting ZYZ
                    dtype          = 'single';
                    angles_refined      = cast(angle,dtype);
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
                    ref_projs = ref_projs(round(size(ref_projs,1)/2)-800+1:round(size(ref_projs,1)/2)+800, ...
                                          round(size(ref_projs,2)/2)-800+1:round(size(ref_projs,2)/2)+800);
                     %% normalize the reference projections
                    ref_projs(ref_projs>0) = ref_projs(ref_projs>0) - min(min(ref_projs(550:950, 400:950)));
                    ref_projs(ref_projs>0) = ref_projs(ref_projs>0) ./max(max(ref_projs(550:950, 400:950)));
                    ref_projs(ref_projs<0) = 0;
                    projs_sm  = projs(round(size(projs,1)/2)-800+1:round(size(projs,1)/2)+800, ...
                                      round(size(projs,2)/2)-800+1:round(size(projs,2)/2)+800,index);
                    mask1 = ref_projs;mask1(mask1>0) = 1;mask1(mask1~=1) = 0;
                    mask2 = projs_sm;mask2(mask2>0) = 1;mask2(mask2~=1) = 0;
                    mask  = mask1.*mask2;
                    ref_projs = ref_projs.*mask;
                    projs_sm  = projs_sm.*mask;
                   
                    
                    %% cross-correlation alignment
                    projs_sm = projs_sm(round(size(projs_sm,1)/2)-797+1:round(size(projs_sm,1)/2)+797, ...
                                        round(size(projs_sm,2)/2)-797+1:round(size(projs_sm,2)/2)+797,:);
                    
                    
                    %% cross correlation to figure out the best match region
                    [colS, rowS] = size(projs_sm);
                    [colB, rowB] = size(ref_projs);
                    corrMatrix = double(zeros(colB - colS + 1, rowB - rowS +1)); 
                    for colPos = 1:colB - colS + 1 
                        for rowPos = 1:rowB - rowS + 1        
                            tempProj = ref_projs(colPos:colPos+colS-1, rowPos:rowPos+rowS-1);
                            corrMatrix(colPos, rowPos) = corr2(tempProj, projs_sm);
                        end
                    end
                    
                    %% find the maximum position
                    [maxNum1, tempR] = max(max(corrMatrix));
                    [maxNum2, tempC] = max(corrMatrix(:,tempR));
                    maxProj =  ref_projs(tempC:tempC+colS-1, tempR:tempR+rowS-1);
                    angle_Matrix(tempx,tempy,tempz) = maxNum2;
                    toc
                    %% update the projection
                    projs(:,:,index) = imtranslate(projs(:,:,index), ...
                        [tempR-round(size(corrMatrix,2)+1)/2,tempC-round(size(corrMatrix,1)+1)/2]);
                end
            end        
        end
        [maxNum0, tempz] = max(max(max(angle_Matrix)));
        [maxNum1, tempy] = max(max(angle_Matrix(:,:,tempz)));
        [maxNum2, tempx] = max(angle_Matrix(:,tempy,tempz));
        angle(1,1) = angle_center(1)+zAxis_angles(tempz);
        angle(1,2) = angle_center(2)+xAxis_angles(tempx);
        angle(1,3) = angle_center(3)+yAxis_angles(tempy);
        save([pwd,'\output_data\Angles_Proj_',num2str(index),'.mat'],'angle');
    
    end
    old_angles = angles;
    for tempi = 1:size(angles,1)
        load([pwd,'\output_data\Angles_Proj_',num2str(tempi),'.mat'])
        angles(tempi,:) = angle;
    end
    if sum(old_angles(:)-angles(:)) == 0
        count = count + 1;
    end
    if count == 6
        break
    end

end
save([pwd,'\output_data\Angles.mat'],'angles');

















