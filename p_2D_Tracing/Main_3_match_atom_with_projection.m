clear;
clc;
%% load data
addpath([pwd,'/src/'])
addpath([pwd,'/input_data/'])
addpath([pwd,'/output_data/'])
projs     = importdata('Projections.mat');
angles    = importdata('Angles.mat');
load('Matched_Flat_Simulation_2D_Tracing.mat')

%% match the 2D Tracing
load('atom_tracing_model_Upper_Layer.mat');clear Slice;atom_up = atom;
load('atom_tracing_model_Lower_Layer.mat');clear Slice;atom_down = atom;
atom = zeros(3,size(atom_up,2)+size(atom_down,2));
atom(1:2,1:size(atom_up,2)) = atom_up;atom(3,1:size(atom_up,2)) = 21;
atom(1:2,size(atom_up,2)+1:size(atom_down,2)+size(atom_up,2)) = atom_down;
atom(3,size(atom_up,2)+1:size(atom_down,2)+size(atom_up,2)) = 57;
atom(3,:) = atom(3,:) - mean(atom(3,:));
for tempi = 1:size(atom_sim,2)
    dif=(atom(:,:)-atom_sim(:,tempi));
    dis=sqrt(sum(dif.^2,1));
    [val,ind] = min(dis);
    if val < 7
        atom_sim(:,tempi) = atom(:,ind);
    end
end

atom_sim = atom_sim.*(0.19/2);

%% projection 2 is the reference
res = 0.19/2;

%% recenter the atom position and crop the projections
atom_sim(1,:) = atom_sim(1,:)-95;
atom_sim(2,:) = atom_sim(2,:)-95;
atom_sim(3,:) = atom_sim(3,:)-mean(atom_sim(3,:));
atom_sim = atom_sim./res;
points = [[-552.2663 -499.0350]
          [-394.9658 -740.4667] 
          [ 227.2067 -713.5666]
          [ 616.1867 -232.0635]
          [ 615.2764  177.9907]
          [ 419.3362  700.5625]
          [  78.5757  701.5758]
          [-210.9211  559.2083]
          [-546.8354  69.9514]];
xv = points(:,1);
yv = points(:,2);
[in, on] = inpolygon(atom_sim(1,:), atom_sim(2,:), xv, yv);
index = in | on;

atomPos = atom_sim(:,index==1);
label   = label_sim(index==1);

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

%% crop the sub region in  projs of each tilting angles
cropProjs = zeros(size( projs));
smallRegion = 18;
[X, Y] = meshgrid(-smallRegion:smallRegion, -smallRegion:smallRegion);                 % Generate the grid
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
    rotScanPos = RotM*flatScanPos(:,:);
    tempCP = zeros([size(cropProjs,1),size(cropProjs,2)]);

    for tempj = 1:size(atomPos,2)

        atomposx = round(rotAtompos(2,tempj));
        atomposy = round(rotAtompos(1,tempj));
        tempCP(atomposy-smallRegion:atomposy+smallRegion,atomposx-smallRegion:atomposx+smallRegion) = ...
            tempCP(atomposy-smallRegion:atomposy+smallRegion,atomposx-smallRegion:atomposx+smallRegion)+cropCircle;
          
    end
    tempCP(tempCP>0) = 1;  
    cropProjs(:,:,tempi) = tempCP.*projs(:,:,tempi);
end

atom = atomPos.*res; 
save([pwd,'/output_data/atom_Refinement_sim_2D_Tracing.mat'],'atom')
save([pwd,'/output_data/label_Refinement_sim_2D_Tracing.mat'],'label')
save([pwd,'/output_data/Projs_Refinement_sim_2D_Tracing.mat'],'cropProjs')
save([pwd,'/output_data/angles_Refinement_sim_2D_Tracing.mat'],'angles')
