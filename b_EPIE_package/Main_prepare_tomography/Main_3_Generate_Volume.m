clear;
clc;
%% load data
addpath([pwd,'\src\'])
addpath([pwd,'\input_data\'])
addpath([pwd,'\output_data\'])
atom_info = importdata('atom_label_Initial_Crop.mat');
para      = importdata('Volume_Parameter.mat');

%% generate the 3D volume model
%% set parameter
atompos   = atom_info.atom;
atomlabel = atom_info.label;
res = 0.19/2;                                                                 % the pixel size unit: angstrom
model = single(atompos);                                                           % atom position coordination
atoms = single(atomlabel);                                                             % atom type
atom_type_num = single(numel(unique(gather(atoms))));
model = single(model./res);

clear label atompos

yMax = single(round(max(model(1,:)) - min(model(1,:)) + 6/(res)));                                                                % the up boundary for y axis
xMax = single(round(max(model(2,:)) - min(model(2,:)) + 6/(res)));                                                                  % the up boundary for x axis
zMax = single(round(max(model(3,:)) - min(model(3,:)) + 6/(res)));   
simulVol = single(zeros(yMax,xMax,zMax));

model(1,:) = single(model(1,:) - (min(model(1,:)) + max(model(1,:)))./2 + round(yMax/2));
model(2,:) = single(model(2,:) - (min(model(2,:)) + max(model(2,:)))./2 + round(xMax/2));
model(3,:) = single(model(3,:) - mean(model(3,:)) + round(zMax/2));        % the z-axis position needed to be centralize 

para = single(reshape(para,[2 atom_type_num]));
peak_value =  single(para(1,:)./para(1,1));   
sigma = single(para(2,:));                                                         % sigma for 3D Gaussian, the first one is for the ordinary atom, the second one is for the marker

% fit the simulation 3D model volumn with atom (using 3D Gaussian as a model of an atom)
for tempj = 1:2
    tempIndex = find((atoms == tempj));
    for tempi = tempIndex
        tic
        [X, Y, Z] = meshgrid(1:xMax, 1:yMax, 1:zMax);
        atomposy = model(1,tempi);        
        atomposx = model(2,tempi);
        atomposz = model(3,tempi);

        tempIy = round(atomposy-30):1:round(atomposy+30);
        tempIx = round(atomposx-30):1:round(atomposx+30);
        tempIz = round(atomposz-30):1:round(atomposz+30);

        tempVol = peak_value(tempj).*exp(-(((X(tempIy,tempIx,tempIz) - atomposx).^2 ...
                                           +(Y(tempIy,tempIx,tempIz) - atomposy).^2 ...
                                           +(Z(tempIy,tempIx,tempIz) - atomposz).^2) / (2*(sigma(tempj)^2))));
        
        simulVol(tempIy,tempIx,tempIz) = simulVol(tempIy,tempIx,tempIz) + tempVol;
        clear tempVol atomposx atomposy atomposz X Y Z
        toc
    end
end
save([pwd,'\output_data\Volume_Atom_Tracing.mat'],'simulVol');
















