function [atom_exp,atom_sim] = atom_align_exp_sim_model_1(atom_exp,atom_sim)
for iteration = 1:2
    %% match the upper and lower region
    atomup = atom_exp(1:2,atom_exp(3,:)<0);
    atomdown = atom_exp(1:2,atom_exp(3,:)>0);
    atomposup = atom_sim(:,atom_sim(3,:)<0);
    atomposdown = atom_sim(:,atom_sim(3,:)>0);
    atomDist = 1.0;
    
    %% re-match the angle
    for iter = 1:2
        %% for the upper region
        if iteration == 1 && iter == 1
            thetaVector = -0.1:0.01:0.1;
        else
            thetaVector = -0.05:0.01:0.05;
        end
        dirVector = zeros(1,size(thetaVector,2));
        count = 1;
        for tempTheta = thetaVector
            tic
            % rotation matrix
            R0 = [cosd(tempTheta) -sind(tempTheta); sind(tempTheta) cosd(tempTheta)];
            tempatomposup = R0*atomposup(1:2,:);
            Mag_shift = 1;
            shift=[0 0]';                                                  % the shift used to do the subpixel shift to align the atom position
            while Mag_shift> 1e-5                                          % if the total difference is larger than 1e-5 (angstrom) then continue to do the alignment
                tempatomposup = tempatomposup  - shift;
                difarr = [];                                               % array for store the difference between two tracing result
                disarr = [];                                               % array for store the distance between two tracing result
                count_arr1 = [];                                           % store the index for the first atom result, which is match with the second atom tracing result
                count_arr2 = [];                                           % store the index for the second atom result, which is match with the first atom tracing result
                for i=1:size(atomup,2)
                    dif=(tempatomposup-atomup(:,i));                       % calculate all difference from the first result with i'th position in the second result (angstrom)
                    dis=sqrt(sum(dif.^2,1));                               % calculate the distance (angstrom)
                    [dis,ind]=min(dis);                                    % obatin the minimum distance and corresponding index
                    if dis <= atomDist                                     % if the minimum distance is smaller than the threshold then store the information
                        difarr=[difarr dif(:,ind)];
                        disarr = [disarr dis];
                    end
                end
                shift=mean(difarr,2);                                      % calculate the mean values of the difference in x-y-z axis, which is also the shift value for the next iteration
                Mag_shift=sum(abs(mean(difarr,2)));                        % calculate the total difference
            end
            dirVector(count) = mean(disarr);
            count = count + 1;
            toc
        end
        [~, index] = min(dirVector);
        thetaUp = thetaVector(index);
        
        %% for the lower region
        dirVector = zeros(1,size(thetaVector,2));
        count = 1;
        for tempTheta = thetaVector
            tic
            % rotation matrix
            R0 = [cosd(tempTheta) -sind(tempTheta); sind(tempTheta) cosd(tempTheta)];
            tempatomposdown = R0*atomposdown(1:2,:);
            Mag_shift = 1;
            shift=[0 0]';                                                  % the shift used to do the subpixel shift to align the atom position
            while Mag_shift> 1e-5                                          % if the total difference is larger than 1e-5 (angstrom) then continue to do the alignment
                tempatomposdown = tempatomposdown  - shift;
                difarr = [];                                               % array for store the difference between two tracing result
                disarr = [];                                               % array for store the distance between two tracing result
                count_arr1 = [];                                           % store the index for the first atom result, which is match with the second atom tracing result
                count_arr2 = [];                                           % store the index for the second atom result, which is match with the first atom tracing result
                for i=1:size(atomdown,2)
                    dif=(tempatomposdown-atomdown(:,i));                   % calculate all difference from the first result with i'th position in the second result (angstrom)
                    dis=sqrt(sum(dif.^2,1));                               % calculate the distance (angstrom)
                    [dis,ind]=min(dis);                                    % obatin the minimum distance and corresponding index
                    if dis <= atomDist                                     % if the minimum distance is smaller than the threshold then store the information
                        difarr=[difarr dif(:,ind)];    
                        disarr = [disarr dis];
                    end
                end
                shift=mean(difarr,2);                                      % calculate the mean values of the difference in x-y-z axis, which is also the shift value for the next iteration
                Mag_shift=sum(abs(mean(difarr,2)));                        % calculate the total difference
            end
            dirVector(count) = mean(disarr);
            count = count + 1;
            toc
        end
        [~, index] = min(dirVector);
        thetaDown  = thetaVector(index);
        if iteration == 1 && iter == 1
            theta      = (thetaDown+thetaUp)./2;
            R0 = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
            atomposdown(1:2,:) = R0*atomposdown(1:2,:);
            atomposup(1:2,:) = R0*atomposup(1:2,:);
        else
            R0 = [cosd(thetaUp) -sind(thetaUp); sind(thetaUp) cosd(thetaUp)];
            atomposup(1:2,:) = R0*atomposup(1:2,:);
            R0 = [cosd(thetaDown) -sind(thetaDown); sind(thetaDown) cosd(thetaDown)];
            atomposdown(1:2,:) = R0*atomposdown(1:2,:);
        end

        %% for upper layer to realign
        Mag_shift = 1;
        shift=[0 0]';                                                      % the shift used to do the subpixel shift to align the atom position
        while Mag_shift> 1e-5                                              % if the total difference is larger than 1e-5 (angstrom) then continue to do the alignment
            atomposup(1:2,:) = atomposup(1:2,:)  - shift;
            difarr = [];                                                   % array for store the difference between two tracing result
            disarr = [];                                                   % array for store the distance between two tracing result
            count_arr1 = [];                                               % store the index for the first atom result, which is match with the second atom tracing result
            count_arr2 = [];                                               % store the index for the second atom result, which is match with the first atom tracing result
            for i=1:size(atomup,2)
                dif=(atomposup(1:2,:)-atomup(:,i));                        % calculate all difference from the first result with i'th position in the second result (angstrom)
                dis=sqrt(sum(dif.^2,1));                                   % calculate the distance (angstrom)
                [dis,ind]=min(dis);                                        % obatin the minimum distance and corresponding index
                if dis <= atomDist                                         % if the minimum distance is smaller than the threshold then store the information
                    difarr=[difarr dif(:,ind)];
                    disarr = [disarr dis];
                end
            end
            shift=mean(difarr,2);                                          % calculate the mean values of the difference in x-y-z axis, which is also the shift value for the next iteration
            Mag_shift=sum(abs(mean(difarr,2)));                            % calculate the total difference
        end
        
        %% for lower layer to realign
        Mag_shift = 1;
        shift=[0 0]';                                                      % the shift used to do the subpixel shift to align the atom position
        while Mag_shift> 1e-5                                              % if the total difference is larger than 1e-5 (angstrom) then continue to do the alignment
            atomposdown(1:2,:) = atomposdown(1:2,:)  - shift;
            difarr = [];                                                   % array for store the difference between two tracing result
            disarr = [];                                                   % array for store the distance between two tracing result
            count_arr1 = [];                                               % store the index for the first atom result, which is match with the second atom tracing result
            count_arr2 = [];                                               % store the index for the second atom result, which is match with the first atom tracing result
            for i=1:size(atomdown,2)
                dif=(atomposdown(1:2,:)-atomdown(:,i));                    % calculate all difference from the first result with i'th position in the second result (angstrom)
                dis=sqrt(sum(dif.^2,1));                                   % calculate the distance (angstrom)
                [dis,ind]=min(dis);                                        % obatin the minimum distance and corresponding index
                if dis <= atomDist                                         % if the minimum distance is smaller than the threshold then store the information
                    difarr=[difarr dif(:,ind)];
                    disarr = [disarr dis];
                end
            end
            shift=mean(difarr,2);                                          % calculate the mean values of the difference in x-y-z axis, which is also the shift value for the next iteration
            Mag_shift=sum(abs(mean(difarr,2)));                            % calculate the total difference
        end
        
        %% re-match the pixel size
        ratioVector = 0.995:0.001:1.005;
        dirVector = zeros(1,size(ratioVector,2));
        count = 1;
        for tempratio = ratioVector
            tempatomposup = tempratio.*atomposup(1:2,:);
            tempatomposdown = tempratio.*atomposdown(1:2,:);
            Mag_shift = 1;
            shift=[0 0]';                                                  % the shift used to do the subpixel shift to align the atom position
            while Mag_shift> 1e-5                                          % if the total difference is larger than 1e-5 (angstrom) then continue to do the alignment
                tempatomposup = tempatomposup  - shift;
                tempatomposdown = tempatomposdown - shift;
                difarr = [];                                               % array for store the difference between two tracing result
                disarr = [];                                               % array for store the distance between two tracing result
                count_arr1 = [];                                           % store the index for the first atom result, which is match with the second atom tracing result
                count_arr2 = [];                                           % store the index for the second atom result, which is match with the first atom tracing result
                for i=1:size(atomup,2)
                    dif=(tempatomposup-atomup(:,i));                       % calculate all difference from the first result with i'th position in the second result (angstrom)
                    dis=sqrt(sum(dif.^2,1));                               % calculate the distance (angstrom)
                    [dis,ind]=min(dis);                                    % obatin the minimum distance and corresponding index
                    if dis <= atomDist                                     % if the minimum distance is smaller than the threshold then store the information
                        difarr=[difarr dif(:,ind)];
                        disarr = [disarr dis];
                    end
                end
                for i=1:size(atomdown,2)
                    dif=(tempatomposdown-atomdown(:,i));                   % calculate all difference from the first result with i'th position in the second result (angstrom)
                    dis=sqrt(sum(dif.^2,1));                               % calculate the distance (angstrom)
                    [dis,ind]=min(dis);                                    % obatin the minimum distance and corresponding index
                    if dis <= atomDist                                     % if the minimum distance is smaller than the threshold then store the information
                        difarr=[difarr dif(:,ind)];
                        disarr = [disarr dis];
                    end
                end
                shift=mean(difarr,2);                                      % calculate the mean values of the difference in x-y-z axis, which is also the shift value for the next iteration
                Mag_shift=sum(abs(mean(difarr,2)));                        % calculate the total difference
            end
            dirVector(count) = mean(disarr);
            count = count + 1;   
        end
    
        [~, index] = min(dirVector);
        ratio = ratioVector(index);
        atomposup(1:2,:) = ratio.*atomposup(1:2,:);
        atomposdown(1:2,:) = ratio.*atomposdown(1:2,:);
        
        Mag_shift = 1;
        sumdisarr = inf;
        shift=[0 0]';                                                      % the shift used to do the subpixel shift to align the atom position
        while Mag_shift> 1e-5                                              % if the total difference is larger than 1e-5 (angstrom) then continue to do the alignment
            atomposup(1:2,:) = atomposup(1:2,:)  - shift;
            difarr = [];                                                   % array for store the difference between two tracing result
            disarr = [];                                                   % array for store the distance between two tracing result
            count_arr1 = [];                                               % store the index for the first atom result, which is match with the second atom tracing result
            count_arr2 = [];                                               % store the index for the second atom result, which is match with the first atom tracing result
            for i=1:size(atomup,2)
                dif=(atomposup(1:2,:)-atomup(:,i));                        % calculate all difference from the first result with i'th position in the second result (angstrom)
                dis=sqrt(sum(dif.^2,1));                                   % calculate the distance (angstrom)
                [dis,ind]=min(dis);                                        % obatin the minimum distance and corresponding index
                if dis <= atomDist                                         % if the minimum distance is smaller than the threshold then store the information
                    difarr=[difarr dif(:,ind)];
                    disarr = [disarr dis];
                end
            end
            if sum(disarr) <= sumdisarr
                sumdisarr = sum(disarr);
            else
                break
            end
            shift=mean(difarr,2);                                          % calculate the mean values of the difference in x-y-z axis, which is also the shift value for the next iteration
            Mag_shift=sum(abs(mean(difarr,2)));                            % calculate the total difference
        end
    
        Mag_shift = 1;
        sumdisarr = inf;
        shift=[0 0]';                                                      % the shift used to do the subpixel shift to align the atom position
        while Mag_shift> 1e-5                                              % if the total difference is larger than 1e-5 (angstrom) then continue to do the alignment
            atomposdown(1:2,:) = atomposdown(1:2,:)  - shift;
            difarr = [];                                                   % array for store the difference between two tracing result
            disarr = [];                                                   % array for store the distance between two tracing result
            count_arr1 = [];                                               % store the index for the first atom result, which is match with the second atom tracing result
            count_arr2 = [];                                               % store the index for the second atom result, which is match with the first atom tracing result
            for i=1:size(atomdown,2)
                dif=(atomposdown(1:2,:)-atomdown(:,i));                    % calculate all difference from the first result with i'th position in the second result (angstrom)
                dis=sqrt(sum(dif.^2,1));                                   % calculate the distance (angstrom)
                [dis,ind]=min(dis);                                        % obatin the minimum distance and corresponding index
                if dis <= atomDist                                         % if the minimum distance is smaller than the threshold then store the information
                    difarr=[difarr dif(:,ind)];
                    disarr = [disarr dis];
                end
            end
            if sum(disarr) <= sumdisarr
                sumdisarr = sum(disarr);
            else
                break
            end
            shift=mean(difarr,2);                                          % calculate the mean values of the difference in x-y-z axis, which is also the shift value for the next iteration
            Mag_shift=sum(abs(mean(difarr,2)));                            % calculate the total difference
        end
    end
    atomposup(3,:) = mean(atom_exp(3,atom_exp(3,:)<0));
    atomposdown(3,:) = mean(atom_exp(3,atom_exp(3,:)>0));
    atom_sim = [atomposup, atomposdown];

    %% crop the same number of simulation atom as experiment atom
    atom_sim_crop = zeros(size(atom_exp));
    atomDist  = 1.0;
    for tempi = 1:size(atom_exp,2)
        dif=(atom_sim-atom_exp(:,tempi));
        dis=sqrt(sum(dif.^2,1));              
        [dis,ind]=min(dis);                   
        if dis <= atomDist
            atom_sim_crop(:,tempi) = atom_sim(:,ind);
        end
    end
    
    %% rotation the experiment to match the flat simulatino model
    while true
        tic
        %% rotation matrix along x-axis
        thetax = -0.05:0.005:0.05;
        meanDifVector = zeros(1,size(thetax,2));
        count = 1;
        for temptethax = thetax
            Rx = [1 0                0                ; 
                  0 cosd(temptethax) -sind(temptethax); 
                  0 sind(temptethax) cosd(temptethax)];
            temp_atom_exp = Rx*atom_exp;
            Mag_shift = 1;
            atomDist  = 1.0 ; % angstrom
            shift=[0 0 0]';
            while Mag_shift> 1e-5                                          % if the total difference is larger than 1e-5 (angstrom) then continue to do the alignment
                temp_atom_exp = temp_atom_exp  - shift;
                difarr = [];                                               % array for store the difference between two tracing result
                disarr = [];                                               % array for store the distance between two tracing result
                for i=1:size(atom_sim_crop,2)
                    dif=(temp_atom_exp-atom_sim_crop(:,i));                % calculate all difference from the first result with i'th position in the second result (angstrom)
                    dis=sqrt(sum(dif.^2,1));                               % calculate the distance (angstrom)
                    [dis,ind]=min(dis);                                    % obatin the minimum distance and corresponding index
                    if dis <= atomDist                                     % if the minimum distance is smaller than the threshold then store the information
                        difarr=[difarr dif(:,ind)];
                        disarr = [disarr dis];
                    end
                end
                shift=mean(difarr,2);                                      % calculate the mean values of the difference in x-y-z axis, which is also the shift value for the next iteration
                Mag_shift=sum(abs(mean(difarr,2)));                        % calculate the total difference
            end
            meanDifVector(count) = mean(disarr);
            count = count + 1;
        end
        [~,ind] = min(meanDifVector);
        rothetax = thetax(ind);
        Rx = [1 0              0              ; 
              0 cosd(rothetax) -sind(rothetax); 
              0 sind(rothetax) cosd(rothetax)];
        atom_exp = Rx*atom_exp;
        %% rotation matrix along y-axis
        thetay = -0.05:0.005:0.05;
        meanDifVector = zeros(1,size(thetay,2));
        count = 1;
        for temptethay = thetay    
            Ry = [cosd(temptethay)  0 sind(temptethay) ; 
                  0                 1 0          ; 
                  -sind(temptethay) 0 cosd(temptethay)];
            temp_atom_exp = Ry*atom_exp;
            Mag_shift = 1;
            atomDist  = 1.0 ;                                              % angstrom
            shift=[0 0 0]';
            while Mag_shift> 1e-5                                          % if the total difference is larger than 1e-5s (angstrom) then continue to do the alignment
                temp_atom_exp = temp_atom_exp  - shift;
                difarr = [];                                               % array for store the difference between two tracing result
                disarr = [];                                               % array for store the distance between two tracing result
                for i=1:size(atom_sim_crop,2)
                    dif=(temp_atom_exp-atom_sim_crop(:,i));                % calculate all difference from the first result with i'th position in the second result (angstrom)
                    dis=sqrt(sum(dif.^2,1));                               % calculate the distance (angstrom)
                    [dis,ind]=min(dis);                                    % obatin the minimum distance and corresponding index
                    if dis <= atomDist                                     % if the minimum distance is smaller than the threshold then store the information
                        difarr=[difarr dif(:,ind)];
                        disarr = [disarr dis];
                    end
                end
                shift=mean(difarr,2);                                      % calculate the mean values of the difference in x-y-z axis, which is also the shift value for the next iteration
                Mag_shift=sum(abs(mean(difarr,2)));                        % calculate the total difference
            end
            meanDifVector(count) = mean(disarr);
            count = count + 1;   
        end
        [~,ind] = min(meanDifVector);
        rothetay = thetay(ind);
        Ry = [cosd(rothetay)  0 sind(rothetay) ; 
              0               1 0          ; 
              -sind(rothetay) 0 cosd(rothetay)];
        atom_exp = Ry*atom_exp;
        if rothetay+rothetax == 0
            toc
            break
        end
        toc
    end
end
end