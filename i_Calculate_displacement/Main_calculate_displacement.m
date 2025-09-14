clear;
clc;

%% load data
addpath([pwd,'/src/'])
addpath([pwd,'/input_data/'])
load('atom_model_sim_exp_matched.mat')

%% load the ,atched atom model data
atom     = [atomup,atomdown];
label    = [labelup,labeldown];
atomgroup= [atomupgroup,atomdowngroup];
atom_sim = [atomup_sim,atomdown_sim];

atomup      = atom(:,atom(3,:)<0);
atomup_sim  = atom_sim(:,atom_sim(3,:)<0); 
atomupgroup = atomgroup(atom(3,:)<0);

atomdown      = atom(:,atom(3,:)>0);
atomdown_sim  = atom_sim(:,atom_sim(3,:)>0); 
atomdowngroup = atomgroup(atom(3,:)>0);

%% check the atom bond distance
atomDist = 190;
atom_1 = [atomup(:,atomupgroup==1),atomdown(:,atomdowngroup==1)];
atom_2 = [atomup(:,atomupgroup==2),atomdown(:,atomdowngroup==2)];
totalAtomdist = zeros(3,size(atom_1,1));
for tempI = 1:size(atom_1,2)
    % count use to count the neighbor atom, which should be equal to after each loop 4 (4 = 1 + 3(3 atom bonds))
    count = 1; 
    for tempJ = 1:size(atom_2,2)
        % calculate the atom bond distance
        deltax = abs(atom_1(1,tempI) - atom_2(1,tempJ));
        deltay = abs(atom_1(2,tempI) - atom_2(2,tempJ));
        deltaz = abs(atom_1(3,tempI) - atom_2(3,tempJ));
        dist =  sqrt(deltax^2+deltay^2+deltaz^2);        
        if dist ~= 0 && dist < atomDist 
            totalAtomdist(count,tempI) = dist;
            count = count + 1;
        end
    end
end
totalDist = [totalAtomdist(1,totalAtomdist(1,:)>0),totalAtomdist(2,totalAtomdist(2,:)>0),totalAtomdist(3,totalAtomdist(3,:)>0)];
bondLength = mean(totalDist);
latticeConstant = bondLength*sqrt(3);
interlayerDistance = abs(mean(atomdown(3,:))-mean(atomup(3,:)));
clear atom_1 atom_2 ax count deltaz deltay deltax dist tempI tempJ totalDist totalAtomdist bondLength

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter setting
convXflag = true;   % if true the x-axis will be convolve
convYflag = true;   % if true the y=axis will be convolve
convZflag = true;   % if true the z-axis will be convolve

sizeNum   = 800;    % size of the Gaussian filter
fraction_xy  = 1;
fraction_z   = 1;
sigmax    = fraction_xy.*latticeConstant;     % sigma for x-axis Gaussian filter
sigmay    = fraction_xy.*latticeConstant;     % sigma for x-axis Gaussian filter
sigmaz    = fraction_z.*interlayerDistance;  % sigma for x-axis Gaussian filter
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% if you want to check the shape of the Gaussian filter please use the code below
%% figure();surf(x2D, y2D, filter(x2D,y2D), 'EdgeColor', 'none');  % Remove grid lines

%% convolve the z-axis
if convZflag == true
    % generate the generalize Gaussian filter
    peak2D     = 1;
    beta2D     = 2;
    [x2D, y2D] = meshgrid(-sizeNum:sizeNum, -sizeNum:sizeNum);
    filter     = double(peak2D.*exp(-(abs(x2D).^beta2D+abs(y2D).^beta2D)./(2*sigmaz^2)));
    filter     = scatteredInterpolant (x2D(:), y2D(:), filter(:),'natural', 'none');
    clear x2D y2D beta2D peak2D sigma latticeConstant 
    
    % convolve the upper layer
    atomNewup = atomup;
    for tempi = 1:size(atomup,2)
        % center current atom coordinates to (0,0) for x and y coordinates 
        tempPos(1,:)  = double(atomup(1,:) - atomup(1,tempi));
        tempPos(2,:)  = double(atomup(2,:) - atomup(2,tempi));
        tempPos(3,:)  = atomup(3,:);
        
        % find the atom inside the region of the Gaussian filter
        tempIndex     = (tempPos(2,:)>-sizeNum) & (tempPos(2,:)<sizeNum) & (tempPos(1,:)>-sizeNum) & (tempPos(1,:)<sizeNum);
        Index         = find(tempIndex == 1);
        sumcoefficent = 0;
        sumZpos       = 0;

        % smooth the z-axis position through weighted average (Similar to Convolution)
        for tempj = 1:size(Index,2)
            coefficent = filter(tempPos(1,Index(tempj)),tempPos(2,Index(tempj)));
            sumcoefficent = sumcoefficent + coefficent;
            sumZpos = sumZpos + tempPos(3,Index(tempj)).*coefficent;
            clear coefficent
        end

        % weighted average and obtain the z coordinates after the convolution
        atomNewup(3,tempi) = sumZpos./sumcoefficent;
        clear Index tempIndex tempPos sumZpos sumZpos
    end
    atomup = atomNewup;
    clear atomNewup

    % convolve the lower layer
    atomNewdown = atomdown;
    for tempi = 1:size(atomdown,2)
        % center current atom coordinates to (0,0) for x and y coordinates 
        tempPos(1,:)  = double(atomdown(1,:) - atomdown(1,tempi));
        tempPos(2,:)  = double(atomdown(2,:) - atomdown(2,tempi));
        tempPos(3,:)  = atomdown(3,:);

        % find the atom inside the region of the Gaussian filter
        tempIndex     = (tempPos(2,:)>-sizeNum) & (tempPos(2,:)<sizeNum) & (tempPos(1,:)>-sizeNum) & (tempPos(1,:)<sizeNum);
        Index         = find(tempIndex == 1);
        sumcoefficent = 0;
        sumZpos       = 0;

        % smooth the z-axis position through weighted average (Similar to Convolution)
        for tempj = 1:size(Index,2)
            coefficent = filter(tempPos(1,Index(tempj)),tempPos(2,Index(tempj)));
            sumcoefficent = sumcoefficent + coefficent;
            sumZpos = sumZpos + tempPos(3,Index(tempj)).*coefficent;
            clear coefficent
        end
        
        % weighted average and obtain the z coordinates after the convolution
        atomNewdown(3,tempi) = sumZpos./sumcoefficent;
        clear Index tempIndex tempPos sumZpos 
    end
    atomdown = atomNewdown;
    clear atomNewdown filter
end

%% convolve the x-axis
if convXflag == true
    % generate the generalize Gaussian filter
    peak2D     = 1;
    beta2D     = 2;
    [x2D, y2D] = meshgrid(-sizeNum:sizeNum, -sizeNum:sizeNum);
    filter     = peak2D.*exp(-(abs(x2D).^beta2D+abs(y2D).^beta2D)./(2*sigmax^2));
    filter     = scatteredInterpolant (x2D(:), y2D(:), filter(:),'natural', 'none');
    clear x2D y2D beta2D peak2D sigma latticeConstant 
    
    % convolve the upper layer
    deltaup_x = atomup_sim(1,:) - atomup(1,:);
    for tempi = 1:size(atomup,2)
        % center current atom coordinates to (0,0) for x and y coordinates 
        tempPos(1,:)  = double(atomup_sim(1,:) - atomup_sim(1,tempi));
        tempPos(2,:)  = double(atomup_sim(2,:) - atomup_sim(2,tempi));
        tempPos(3,:)  = deltaup_x;

        % find the atom inside the region of the Gaussian filter
        tempIndex     = (tempPos(2,:)>-sizeNum) & (tempPos(2,:)<sizeNum) & (tempPos(1,:)>-sizeNum) & (tempPos(1,:)<sizeNum);
        Index         = find(tempIndex == 1);
        sumcoefficent = 0;
        sumXpos       = 0;

        % smooth the x-axis position through weighted average (Similar to Convolution)
        for tempj = 1:size(Index,2)
            coefficent = filter(tempPos(1,Index(tempj)),tempPos(2,Index(tempj)));
            sumcoefficent = sumcoefficent + coefficent;
            sumXpos = sumXpos + tempPos(3,Index(tempj)).*coefficent;
            clear coefficent
        end
        
        % weighted average and obtain the x coordinates after the convolution
        atomup(1,tempi) = atomup_sim(1,tempi) - sumXpos./sumcoefficent;
        clear Index tempIndex tempPos sumXpos 
    end
    clear  deltaup_x

    % convolve the lower layer
    deltadown_x = atomdown_sim(1,:) - atomdown(1,:);
    for tempi = 1:size(atomdown,2)
        % center current atom coordinates to (0,0) for x and y coordinates 
        tempPos(1,:)  = double(atomdown_sim(1,:) - atomdown_sim(1,tempi));
        tempPos(2,:)  = double(atomdown_sim(2,:) - atomdown_sim(2,tempi));
        tempPos(3,:)  = deltadown_x;

        % find the atom inside the region of the Gaussian filter
        tempIndex     = (tempPos(2,:)>-sizeNum) & (tempPos(2,:)<sizeNum) & (tempPos(1,:)>-sizeNum) & (tempPos(1,:)<sizeNum);
        Index         = find(tempIndex == 1);
        sumcoefficent = 0;
        sumXpos       = 0;

        % smooth the x-axis position through weighted average (Similar to Convolution)
        for tempj = 1:size(Index,2)
            coefficent = filter(tempPos(1,Index(tempj)),tempPos(2,Index(tempj)));
            sumcoefficent = sumcoefficent + coefficent;
            sumXpos = sumXpos + tempPos(3,Index(tempj)).*coefficent;
            clear coefficent
        end
        
        % weighted average and obtain the x coordinates after the convolution
        atomdown(1,tempi) = atomdown_sim(1,tempi) - sumXpos./sumcoefficent;
        clear Index tempIndex tempPos sumXpos 
    end
    clear  deltadown_x filter
end

%% convolve the y-axis
if convXflag == true
    % generate the generalize Gaussian filter
    peak2D     = 1;
    beta2D     = 2;
    [x2D, y2D] = meshgrid(-sizeNum:sizeNum, -sizeNum:sizeNum);
    filter     = peak2D.*exp(-(abs(x2D).^beta2D+abs(y2D).^beta2D)./(2*sigmay^2));
    filter     = scatteredInterpolant (x2D(:), y2D(:), filter(:),'natural', 'none');
    clear x2D y2D beta2D peak2D sigma latticeConstant 
    
    % convolve the upper layer
    deltaup_y = atomup_sim(2,:) - atomup(2,:);
    for tempi = 1:size(atomup,2)
        % center current atom coordinates to (0,0) for x and y coordinates 
        tempPos(1,:)  = double(atomup_sim(1,:) - atomup_sim(1,tempi));
        tempPos(2,:)  = double(atomup_sim(2,:) - atomup_sim(2,tempi));
        tempPos(3,:)  = deltaup_y;

        % find the atom inside the region of the Gaussian filter
        tempIndex     = (tempPos(2,:)>-sizeNum) & (tempPos(2,:)<sizeNum) & (tempPos(1,:)>-sizeNum) & (tempPos(1,:)<sizeNum);
        Index         = find(tempIndex == 1);
        sumcoefficent = 0;
        sumYpos       = 0;

        % smooth the y-axis position through weighted average (Similar to Convolution)
        for tempj = 1:size(Index,2)
            coefficent = filter(tempPos(1,Index(tempj)),tempPos(2,Index(tempj)));
            sumcoefficent = sumcoefficent + coefficent;
            sumYpos = sumYpos + tempPos(3,Index(tempj)).*coefficent;
            clear coefficent
        end
        
        % weighted average and obtain the x coordinates after the convolution
        atomup(2,tempi) = atomup_sim(2,tempi) - sumYpos./sumcoefficent;
        clear Index tempIndex tempPos sumYpos 
    end
    clear deltaup_y
    
    % convolve the lower layer
    deltadown_y = atomdown_sim(2,:) - atomdown(2,:);
    for tempi = 1:size(atomdown,2)
        % center current atom coordinates to (0,0) for x and y coordinates 
        tempPos(1,:)  = double(atomdown_sim(1,:) - atomdown_sim(1,tempi));
        tempPos(2,:)  = double(atomdown_sim(2,:) - atomdown_sim(2,tempi));
        tempPos(3,:)  = deltadown_y;

        % find the atom inside the region of the Gaussian filter
        tempIndex     = (tempPos(2,:)>-sizeNum) & (tempPos(2,:)<sizeNum) & (tempPos(1,:)>-sizeNum) & (tempPos(1,:)<sizeNum);
        Index         = find(tempIndex == 1);
        sumcoefficent = 0;
        sumYpos       = 0;

        % smooth the y-axis position through weighted average (Similar to Convolution)
        for tempj = 1:size(Index,2)
            coefficent = filter(tempPos(1,Index(tempj)),tempPos(2,Index(tempj)));
            sumcoefficent = sumcoefficent + coefficent;
            sumYpos = sumYpos + tempPos(3,Index(tempj)).*coefficent;
            clear coefficent
        end
        
        % weighted average and obtain the x coordinates after the convolution
        atomdown(2,tempi) = atomdown_sim(2,tempi) - sumYpos./sumcoefficent;
        clear Index tempIndex tempPos sumYpos 
    end
    clear deltadown_y filter
end
% realign the z-axis again
atomup_sim(3,:)   = mean(atomup(3,:));
atomdown_sim(3,:) = mean(atomdown(3,:));

%% calculate the displacement and display the vector in upper layer
displacementUP = atomup - atomup_sim;
figure(1)
quiver3(atomup_sim(1,:), atomup_sim(2,:), atomup_sim(3,:), displacementUP(1,:).*50, displacementUP(2,:).*50, displacementUP(3,:).*50 ,'Color', 'b', 'LineWidth', 1, 'MaxHeadSize', 4);
axis equal;
view(90,30);
hold off;
save([pwd,'/output_data/atom_model_upper_layer_displacement_1_convolution.mat'],'atomup_sim','displacementUP')

%% calculate the displacement and display the vector in lower layer
displacementDOWN = atomdown - atomdown_sim;
figure(2)
quiver3(atomdown_sim(1,:), atomdown_sim(2,:), atomdown_sim(3,:), displacementDOWN(1,:).*50, displacementDOWN(2,:).*50, displacementDOWN(3,:).*50 ,'Color', 'r', 'LineWidth', 1, 'MaxHeadSize', 4);
axis equal;
view(90,30);
hold off;
save([pwd,'/output_data/atom_model_lower_layer_displacement_1_Convolution.mat'],'atomdown_sim','displacementDOWN')
