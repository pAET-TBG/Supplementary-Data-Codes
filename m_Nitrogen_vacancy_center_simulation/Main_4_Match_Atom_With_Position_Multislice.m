clear;
clc;
%% load data
addpath([pwd,'/src/'])
addpath([pwd,'/input_data/'])
addpath([pwd,'/output_data/'])
projs     = importdata('Projs_NVCenter_Multislice.mat');
angles    = importdata('Angles_NVcenter.mat');
atompos   = importdata('atom_manual_tracing_NVC_Multislice_With_Probe_Total.mat');
label     = ones(1,size(atompos,2));

%% projection 2 is the reference
res = 0.2129;

%% recenter the atom position and crop the projections
atompos(1,:) = atompos(1,:)-mean(atompos(1,:))-0.05;
atompos(2,:) = atompos(2,:)-mean(atompos(2,:))+0.1;
atompos(3,:) = atompos(3,:)-mean(atompos(3,:));
atomPos = atompos./res;

%% calculate the scan position
halfColnum = 0;                                                            % the half number of the column scan position
halfRownum = 0;                                                            % the half number of the row scan position
% flatScanPos is a 3-by-N matrix, (1,:) is the y-position, (2,:) is the x-position, (3,:) is the z-position
flatScanPos = zeros(3,(halfColnum*2+1)*(halfRownum*2+1));                  % initialize the scan position in the flat plane without tilting (matrix size: (3, N))
flatScanPos(3,:) = 0;                                                      % the z-axis of the scan position in the flat plane without tilting is zero
[tempx, tempy] = meshgrid(-halfRownum:halfRownum, -halfColnum:halfColnum);
flatScanPos(1,:) = reshape(tempy, 1, size(flatScanPos(1,:),2));
flatScanPos(2,:) = reshape(tempx, 1, size(flatScanPos(2,:),2));
flatScanPos(1,:) = flatScanPos(1,:).*(size(projs,1))./(halfColnum*2+1);
flatScanPos(2,:) = flatScanPos(2,:).*(size(projs,2))./(halfRownum*2+1);
clear tempx tempy

%% calculate the scan position for each tilting angles
scanPos = zeros(2,(halfColnum*2+1)*(halfRownum*2+1),size(angles,1));       % initialize the scan position for all tilting angles (matrix size: (2, N, number of tilting angles)) 
for tempi = 1:size(angles,1)
    % calculate the rotation matrix
    R1 = MatrixQuaternionRot([0 0 1],angles(tempi,1));  
    R2 = MatrixQuaternionRot([0 1 0],angles(tempi,2));
    R3 = MatrixQuaternionRot([1 0 0],angles(tempi,3));  
    RotM   = (R1*R2*R3)';
    rotScanPos = RotM*flatScanPos;                                         % the rotated scan position
    rotScanPos(1,:) = rotScanPos(1,:) + (size( projs,1)+1)/2;
    rotScanPos(2,:) = rotScanPos(2,:) + (size( projs,2)+1)/2;    
    scanPos(1,:,tempi) = rotScanPos(1,:); 
    scanPos(2,:,tempi) = rotScanPos(2,:);  
    clear rotScanPos  RotM 
end

%% check position for each  projs 
%% crop the sub region in  projs of each tilting angles 
cropProjs = zeros(size( projs));
smallRegion = 18;
[X, Y] = meshgrid(-smallRegion:smallRegion, -smallRegion:smallRegion);     % Generate the grid
cropCircle = sqrt(X.^2 + Y.^2);
cropCircle(cropCircle<=smallRegion) = 1;cropCircle(cropCircle>smallRegion) = 0;

% the code below is used to check whether the atom position is correct
for tempi = 1:size(projs,3) 
    Rotangle = angles(tempi,:);                                            % corresponding rotation angle
    % calculate the rotation matrix
    RotM = (MatrixQuaternionRot([0;0;1],Rotangle(1))*MatrixQuaternionRot([0;1;0],Rotangle(2))*MatrixQuaternionRot([1;0;0],Rotangle(3)))';
    rotAtompos = RotM*atomPos;                                             % the rotated scan position
    rotAtompos(1,:) = rotAtompos(1,:) + (size(projs(:,:,tempi),1)+1)/2;
    rotAtompos(2,:) = rotAtompos(2,:) + (size(projs(:,:,tempi),2)+1)/2;
end

atom = atomPos.*res; 
save([pwd,'/output_data/atom_Tracing_NVC_Multislice_With_Probe_Total.mat'],'atom')
save([pwd,'/output_data/label_Tracing_NVC_Multislice_With_Probe_Total.mat'],'label')
save([pwd,'/output_data/Projs_Tracing_NVC_Multislice_With_Probe_Total.mat'],'projs')
save([pwd,'/output_data/angles_Tracing_NVC_Multislice_With_Probe_Total.mat'],'angles')

