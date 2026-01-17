function [atomup,atomup_sim,atomdown,atomdown_sim] = atom_align_exp_sim_model_2(atomup,atomup_sim,atomdown,atomdown_sim)
atom_exp  = [atomup,atomdown];
atom_sim  = [atomup_sim,atomdown_sim];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% THE UNIT OF THE COORDINATES IS (pm) %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% center the coordinates for rotation 
atom_exp(1,:) = atom_exp(1,:) - mean(atom_sim(1,:));
atom_exp(2,:) = atom_exp(2,:) - mean(atom_sim(2,:));
atom_exp(3,:) = atom_exp(3,:) - mean(atom_sim(3,:));
atom_sim(1,:) = atom_sim(1,:) - mean(atom_sim(1,:));
atom_sim(2,:) = atom_sim(2,:) - mean(atom_sim(2,:));
atom_sim(3,:) = atom_sim(3,:) - mean(atom_sim(3,:));

for iteration = 1:2
    %% match the upper and lower region separately
    atomup = atom_exp(1:2,atom_exp(3,:)<0);
    atomdown = atom_exp(1:2,atom_exp(3,:)>0);
    atomposup = atom_sim(1:2,atom_sim(3,:)<0);
    atomposdown = atom_sim(1:2,atom_sim(3,:)>0);         
    
    %% re-match the angle (this rotation is along z axis)
    for iter = 1:3
        %% for the upper region
        if iteration == 1 
            thetaVector = -0.3:0.001:0.3;
        else
            thetaVector = -0.05:0.0001:0.05;
        end
        disVector = zeros(1,size(thetaVector,2));                          % the parameter to the average distance value between the flat model and experiment coordinates
        indexCount = 1;                                                    % count the index for the disVector
        for tempTheta = thetaVector
            tic
            % rotation matrix
            R0 = [cosd(tempTheta) -sind(tempTheta); sind(tempTheta) cosd(tempTheta)];
            tempatomposup = R0*atomposup(1:2,:);                            
            Mag_shift = 1;
            shift=[0 0]';                                                  % the shift used to do the subpixel shift to align the atom position
            while Mag_shift> 1e-6                                          % if the total difference is larger than 1e-6 (angstrom) then continue to do the alignment
                tempatomposup = tempatomposup - shift;                     % shift the flat model
                difarr    = tempatomposup - atomup;                        % calculate the displacement vector
                disarr    = sqrt(sum(difarr.^2,1));                        % calculate the total displacement distance length
                shift     = mean(difarr,2);                                % calculate the mean values of the difference in x-y-z axis, which is also the shift value for the next iteration
                Mag_shift =sum(abs(mean(difarr,2)));                       % calculate the total vector difference
            end
            disVector(indexCount) = mean(disarr);                          % store the average distance value between the flat model and experiment coordinates
            indexCount = indexCount + 1;
            toc
        end
        [~, index] = min(disVector);                                       % choose the minimized value
        thetaUp = thetaVector(index);                                      % set the rotation angle
        
        %% for the lower region
        disVector = zeros(1,size(thetaVector,2));                          % the parameter to the average distance value between the flat model and experiment coordinates
        indexCount = 1;                                                    % count the index for the disVector
        for tempTheta = thetaVector
            tic
            % rotation matrix
            R0 = [cosd(tempTheta) -sind(tempTheta); sind(tempTheta) cosd(tempTheta)];
            tempatomposdown = R0*atomposdown(1:2,:);
            Mag_shift = 1;
            shift=[0 0]';                                                  % the shift used to do the subpixel shift to align the atom position
            while Mag_shift> 1e-6                                          % if the total difference is larger than 1e-6 (angstrom) then continue to do the alignment
                tempatomposdown = tempatomposdown - shift;                 % shift the flat model
                difarr    = tempatomposdown - atomdown;                    % calculate the displacement vector
                disarr    = sqrt(sum(difarr.^2,1));                        % calculate the total displacement distance length
                shift     = mean(difarr,2);                                % calculate the mean values of the difference in x-y-z axis, which is also the shift value for the next iteration
                Mag_shift =sum(abs(mean(difarr,2)));                       % calculate the total vector difference
            end
            disVector(indexCount) = mean(disarr);
            indexCount = indexCount + 1;
            toc
        end
        [~, index] = min(disVector);                                       % choose the minimized value
        thetaDown  = thetaVector(index);                                   % set the rotation angle

        % rotate the two layer
        R0 = [cosd(thetaUp) -sind(thetaUp); sind(thetaUp) cosd(thetaUp)];
        atomposup(1:2,:) = R0*atomposup(1:2,:);
        R0 = [cosd(thetaDown) -sind(thetaDown); sind(thetaDown) cosd(thetaDown)];
        atomposdown(1:2,:) = R0*atomposdown(1:2,:);

        %% for upper layer to realign
        Mag_shift = 1;
        shift=[0 0]';                                                      % the shift used to do the subpixel shift to align the atom position
        while Mag_shift> 1e-6                                              % if the total difference is larger than 1e-6 (angstrom) then continue to do the alignment
            atomposup(1:2,:) = atomposup(1:2,:)  - shift;                  % shift the flat model
            difarr    = atomposup(1:2,:) - atomup;                         % calculate the displacement vector
            disarr    = sqrt(sum(difarr.^2,1));                            % calculate the total displacement distance length
            shift     = mean(difarr,2);                                    % calculate the mean values of the difference in x-y-z axis, which is also the shift value for the next iteration
            Mag_shift =sum(abs(mean(difarr,2)));                           % calculate the total vector difference
        end
        
        %% for lower layer to realign
        Mag_shift = 1;
        shift=[0 0]';                                                      % the shift used to do the subpixel shift to align the atom position
        while Mag_shift> 1e-6                                              % if the total difference is larger than 1e-6 (angstrom) then continue to do the alignment
            atomposdown(1:2,:) = atomposdown(1:2,:)  - shift;              % shift the flat model
            difarr    = atomposdown(1:2,:) - atomdown;                     % calculate the displacement vector
            disarr    = sqrt(sum(difarr.^2,1));                            % calculate the total displacement distance length
            shift     = mean(difarr,2);                                    % calculate the mean values of the difference in x-y-z axis, which is also the shift value for the next iteration
            Mag_shift =sum(abs(mean(difarr,2)));                           % calculate the total vector difference
        end
        
        %% re-match the pixel size
        ratioVector = 0.995:0.0001:1.005;                                  
        disVector = zeros(1,size(ratioVector,2));                          % the parameter to the average distance value between the flat model and experiment coordinates
        indexCount = 1;                                                    % count the index for the disVector
        for tempratio = ratioVector
            % change the size of the upper and lower layer simultaneously 
            tempatomposup = tempratio.*atomposup(1:2,:);                   
            tempatomposdown = tempratio.*atomposdown(1:2,:);
            Mag_shift = 1;
            shift=[0 0]';                                                  % the shift used to do the subpixel shift to align the atom position
            while Mag_shift> 1e-6                                          % if the total difference is larger than 1e-6 (angstrom) then continue to do the alignment
                tempatomposup = tempatomposup  - shift;                    % shift the flat model
                tempatomposdown = tempatomposdown - shift;                 % shift the flat model
                % calculate the displacement vector
                difarr    = [tempatomposup,tempatomposdown] - [atomup,atomdown];
                disarr    = sqrt(sum(difarr.^2,1));                        % calculate the total displacement distance length
                shift     = mean(difarr,2);                                % calculate the mean values of the difference in x-y-z axis, which is also the shift value for the next iteration
                Mag_shift =sum(abs(mean(difarr,2)));                       % calculate the total vector difference        
            end
            disVector(indexCount) = mean(disarr);                          % choose the minimized value
            indexCount = indexCount + 1;                                   % set the rotation angle   
        end
    
        [~, index] = min(disVector);                                       % choose the minimized value
        ratio = ratioVector(index);                                        % set the ratio
        
        % change the size of the two layer
        atomposup(1:2,:) = ratio.*atomposup(1:2,:);
        atomposdown(1:2,:) = ratio.*atomposdown(1:2,:);

        Mag_shift = 1;
        shift=[0 0]';                                                      % the shift used to do the subpixel shift to align the atom position
        while Mag_shift> 1e-6                                              % if the total difference is larger than 1e-6 (angstrom) then continue to do the alignment
            atomposup(1:2,:) = atomposup(1:2,:)  - shift;                  % shift the flat model
            difarr    = atomposup(1:2,:) - atomup;                         % calculate the displacement vector
            disarr    = sqrt(sum(difarr.^2,1));                            % calculate the total displacement distance length
            shift     = mean(difarr,2);                                    % calculate the mean values of the difference in x-y-z axis, which is also the shift value for the next iteration
            Mag_shift =sum(abs(mean(difarr,2)));                           % calculate the total vector difference   
        end
    
        Mag_shift = 1;
        shift=[0 0]';                                                      % the shift used to do the subpixel shift to align the atom position
        while Mag_shift> 1e-6                                              % if the total difference is larger than 1e-6 (angstrom) then continue to do the alignment
            atomposdown(1:2,:) = atomposdown(1:2,:)  - shift;              % shift the flat model
            difarr    = atomposdown(1:2,:) - atomdown;                     % calculate the displacement vector
            disarr    = sqrt(sum(difarr.^2,1));                            % calculate the total displacement distance length
            shift     = mean(difarr,2);                                    % calculate the mean values of the difference in x-y-z axis, which is also the shift value for the next iteration
            Mag_shift =sum(abs(mean(difarr,2)));                           % calculate the total vector difference   
        end
    end

    % set the z coordinates for the flat model as the average value of the experiment z coordinates
    atomposup(3,:) = mean(atom_exp(3,atom_exp(3,:)<0));
    atomposdown(3,:) = mean(atom_exp(3,atom_exp(3,:)>0));
    % store the new fiting result
    atom_sim = [atomposup, atomposdown];
 
    %% rotation the experiment to match the flat simulatino model (this rotation is along x and y axis)
    while true
        tic
        %% rotation matrix along x-axis
        thetax = -0.05:0.0001:0.05;
        disVector = zeros(1,size(thetax,2));                               % the parameter to the average distance value between the flat model and experiment coordinates
        indexCount = 1;                                                    % count the index for the disVector
        for temptethax = thetax
            % rotation matrix
            Rx = [1 0                0                ; 
                  0 cosd(temptethax) -sind(temptethax); 
                  0 sind(temptethax) cosd(temptethax)];
            temp_atom_exp = Rx*atom_exp;
            Mag_shift = 1;
            shift=[0 0 0]';                                                % the shift used to do the subpixel shift to align the atom position
            while Mag_shift> 1e-5                                          % if the total difference is larger than 1e-6 (angstrom) then continue to do the alignment
                temp_atom_exp = temp_atom_exp  - shift;                    % shift the flat model
                difarr    = temp_atom_exp - atom_sim;                      % calculate the displacement vector
                disarr    = sqrt(sum(difarr.^2,1));                        % calculate the total displacement distance length
                shift     = mean(difarr,2);                                % calculate the mean values of the difference in x-y-z axis, which is also the shift value for the next iteration
                Mag_shift =sum(abs(mean(difarr,2)));                       % calculate the total vector difference
            end
            disVector(indexCount) = mean(disarr);
            indexCount = indexCount + 1;
        end
        [~,ind] = min(disVector);                                          % choose the minimized value
        % rotate the two layer
        rothetax = thetax(ind);
        Rx = [1 0              0              ; 
              0 cosd(rothetax) -sind(rothetax); 
              0 sind(rothetax) cosd(rothetax)];
        atom_exp = Rx*atom_exp;
        
        %% rotation matrix along y-axis
        thetay = -0.05:0.0001:0.05;
        disVector = zeros(1,size(thetay,2));                               % the parameter to the average distance value between the flat model and experiment coordinates
        indexCount = 1;                                                    % count the index for the disVector
        for temptethay = thetay    
            % rotation matrix
            Ry = [cosd(temptethay)  0 sind(temptethay) ; 
                  0                 1 0          ; 
                  -sind(temptethay) 0 cosd(temptethay)];
            temp_atom_exp = Ry*atom_exp;
            Mag_shift = 1;
            shift=[0 0 0]';                                                % the shift used to do the subpixel shift to align the atom position
            while Mag_shift> 1e-5                                          % if the total difference is larger than 1e-6 (angstrom) then continue to do the alignment
                temp_atom_exp = temp_atom_exp  - shift;                    % shift the flat model
                difarr    = temp_atom_exp - atom_sim;                      % calculate the displacement vector
                disarr    = sqrt(sum(difarr.^2,1));                        % calculate the total displacement distance length
                shift     = mean(difarr,2);                                % calculate the mean values of the difference in x-y-z axis, which is also the shift value for the next iteration
                Mag_shift =sum(abs(mean(difarr,2)));                       % calculate the total difference
            end
            disVector(indexCount) = mean(disarr);
            indexCount = indexCount + 1;   
        end
        [~,ind] = min(disVector);                                          % choose the minimized value
        % rotate the two layer
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
atomup      = atom_exp(:,atom_exp(3,:)<0);
atomup_sim  = atom_sim(:,atom_sim(3,:)<0); 

atomdown      = atom_exp(:,atom_exp(3,:)>0);
atomdown_sim  = atom_sim(:,atom_sim(3,:)>0); 

end