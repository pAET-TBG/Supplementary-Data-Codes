clear;
clc;

%% load data
addpath([pwd,'\src\'])
addpath([pwd,'\input_data\'])
addpath([pwd,'\output_data\'])

for idx = 1:2
    if idx == 1
        load('atom_tracing_model_original_matched.mat')
    elseif idx == 2
        load('atom_tracing_model_refinement_total_matched.mat')
    end

    %% load the ,atched atom model data
    atom     = [atomup,atomdown];
    atom_sim = [atomup_sim,atomdown_sim];
    atom_match = zeros(size(atom));
    atomDist = 60;
    Mag_shift=1;
    shift=[0 0 0]';   % the shift used to do the subpixel shift to align the atom position
    while Mag_shift>1e-5 % if the total difference is larger than 1e-6 (angstrom) then continue to do the alignment
        atom =  atom - shift; % substruct the difference (i.e. do the shift to align the atom)
        difarr = [];    % array for store the difference between two tracing result
        for i=1:size(atom_sim,2)
            dif=(atom-atom_sim(:,i));         % calculate all difference from the first result with i'th position in the second result (angstrom)
            dis=sqrt(sum(dif.^2,1));                  % calculate the distance (angstrom)
            [dis,ind]=min(dis);                       % obatin the minimum distance and corresponding index
            if dis <= atomDist                        % if the minimum distance is smaller than the threshold then store the information
                difarr=[difarr dif(:,ind)];
                atom_match(:,i) = atom(:,ind);
            end
        end
        shift=mean(difarr,2);                  % calculate the mean values of the difference in x-y-z axis, which is also the shift value for the next iteration
        Mag_shift=sum(abs(mean(difarr,2)));          % calculate the total difference
    end
    clear Mag_shift shift difarr shift dif dis i ind atomDist  
    atom = atom_match;
    
    atomup      = atom(:,atom(3,:)<0);
    atomup_sim  = atom_sim(:,atom_sim(3,:)<0); 
    
    atomdown      = atom(:,atom(3,:)>0);
    atomdown_sim  = atom_sim(:,atom_sim(3,:)>0);  

    %% group the atom
    atomupgroup = zeros(1,size(atomup,2));
    atomupgroup(1) = 1;inside_1 = 1;inside_2 = [];
    while sum(atomupgroup==0) > 0
        for tempi = inside_1
            dif=(atomup-atomup(:,tempi));                    % calculate all difference from the first result with i'th position in the second result (angstrom)
            dis=sqrt(sum(dif.^2,1));    
            inside_2 = [inside_2 find(dis>0&dis<190&atomupgroup==0)];
            atomupgroup(dis>0&dis<190) = 2;
            inside_1(1) = [];
        end    
        for tempi = inside_2
            dif=(atomup-atomup(:,tempi));                    % calculate all difference from the first result with i'th position in the second result (angstrom)
            dis=sqrt(sum(dif.^2,1));    
            inside_1 = [inside_1 find(dis>0&dis<190&atomupgroup==0)];
            atomupgroup(dis>0&dis<190) = 1;
            inside_2(1) = [];
        end
    end
    atomdowngroup = zeros(1,size(atomdown,2));
    atomdowngroup(1) = 1;inside_1 = 1;inside_2 = [];
    while sum(atomdowngroup==0) > 0
        for tempi = inside_1
            dif=(atomdown-atomdown(:,tempi));                    % calculate all difference from the first result with i'th position in the second result (angstrom)
            dis=sqrt(sum(dif.^2,1));    
            inside_2 = [inside_2 find(dis>0&dis<190&atomdowngroup==0)];
            atomdowngroup(dis>0&dis<190) = 2;
            inside_1(1) = [];
        end    
        for tempi = inside_2
            dif=(atomdown-atomdown(:,tempi));                    % calculate all difference from the first result with i'th position in the second result (angstrom)
            dis=sqrt(sum(dif.^2,1));    
            inside_1 = [inside_1 find(dis>0&dis<190&atomdowngroup==0)];
            atomdowngroup(dis>0&dis<190) = 1;
            inside_2(1) = [];
        end
    end
    atomDist = 190;
    atom_1 = [atomup(:,atomupgroup==1),atomdown(:,atomdowngroup==1)];
    atom_2 = [atomup(:,atomupgroup==2),atomdown(:,atomdowngroup==2)];
    totalAtomdist = zeros(3,size(atom_1,2));
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
    figure(12)
    if idx == 1
        subplot(1,2,1)
        quiver3(atomup_sim(1,:), atomup_sim(2,:), atomup_sim(3,:), displacementUP(1,:).*30, displacementUP(2,:).*30, displacementUP(3,:).*30,0 ,'Color', 'b', 'LineWidth', 1, 'MaxHeadSize', 4,'AutoScale', 'off');hold on
        view(0,90);
        xlim([-3500 3500])
        ylim([-5500 5500])
        zlim([-450 -000])
        axis equal tight off;
        title('Ground Truth')
        hold off;
    elseif idx == 2
        subplot(1,2,2)
        quiver3(atomup_sim(1,:), atomup_sim(2,:), atomup_sim(3,:), displacementUP(1,:).*30, displacementUP(2,:).*30, displacementUP(3,:).*30,0 ,'Color', 'b', 'LineWidth', 1, 'MaxHeadSize', 4,'AutoScale', 'off');hold on
        view(0,90);
        xlim([-3500 3500])
        ylim([-5500 5500])
        zlim([-450 -000])
        axis equal tight off;
        title('Traced Result')
        hold off;
    end
    
    %% calculate the displacement and display the vector in lower layer
    displacementDOWN = atomdown - atomdown_sim;
    figure(22)
    if idx == 1
        subplot(1,2,1)
        quiver3(atomdown_sim(1,:), atomdown_sim(2,:), atomdown_sim(3,:), displacementDOWN(1,:).*30, displacementDOWN(2,:).*30, displacementDOWN(3,:).*30,0 ,'Color', 'r', 'LineWidth', 1, 'MaxHeadSize', 4,'AutoScale', 'off');
        view(0,90);
        xlim([-3500 3500])
        ylim([-5500 5500])
        zlim([000 450])
        axis equal tight off;
        title('Ground Truth')
        hold off;
    elseif idx == 2
        subplot(1,2,2)
        quiver3(atomdown_sim(1,:), atomdown_sim(2,:), atomdown_sim(3,:), displacementDOWN(1,:).*30, displacementDOWN(2,:).*30, displacementDOWN(3,:).*30,0 ,'Color', 'r', 'LineWidth', 1, 'MaxHeadSize', 4,'AutoScale', 'off');
        view(0,90);
        xlim([-3500 3500])
        ylim([-5500 5500])
        axis equal tight off;
        title('Traced Result')
        zlim([000 450])
        hold off;
    end

end
