clear;clc;

%% load the data
addpath([pwd,'/src/'])
addpath([pwd,'/input_data/'])
addpath([pwd,'/output_data/'])
atom_info = importdata('atom_model_upper_layer_displacement_1_Convolution.mat');
atom_sim_pos  = atom_info.atomup_sim;
atom_displace = atom_info.displacementUP;
atom_info = importdata('atom_model_sim_exp_matched.mat');
atom_group = atom_info.atomupgroup;
clear atom_info

%% rotate the coordinates
atom_head = atom_sim_pos + atom_displace;
% Define the angle of rotation in degrees
angle_degrees = 180;
% Rotation matrix for x-axis
Rx = [1 0 0;
      0 cosd(angle_degrees) -sind(angle_degrees);
      0 sind(angle_degrees) cosd(angle_degrees)];
atom_sim_pos  = Rx*atom_sim_pos;
atom_head     = Rx*atom_head;
atom_displace = atom_head - atom_sim_pos;
atom_displace_old = atom_displace;
clear atom_head angle_degrees 

%% normalize the displacement vector
displaceNorm = sqrt(sum(atom_displace.^2,1));
atom_displace(1,:) = atom_displace(1,:)./displaceNorm;
atom_displace(2,:) = atom_displace(2,:)./displaceNorm;
atom_displace(3,:) = atom_displace(3,:)./displaceNorm;
clear displaceNorm

%% caculation the angle and group the vector
theta_arr = zeros(1,size(atom_sim_pos,2));
for tempi = 1:size(theta_arr,2)
    lenH = sqrt(atom_displace(1,tempi).^2 + atom_displace(2,tempi).^2);
    lenV = atom_displace(3,tempi);
    theta_arr(tempi) = atand(lenV/lenH);
    clear lenV lenH
end
clear tempi
index_group = theta_arr < 10;
atom_displace_old = atom_displace_old(:,index_group);
atom_displace = atom_displace(:,index_group);
atom_sim_pos = atom_sim_pos(:,index_group);
atom_group = atom_group(:,index_group);
clear index_group


%% form the triangle region sample
sample_points = [  875.5001 -5.2528e+03 0                                    % index 1 2 3 is a triangle region
                 1.0048e+03 -5.3164e+03 0                                    % index 1 3 4 is a triangle region
                 1.1245e+03 -5.2362e+03 0                                    % index 1 4 5 is a triangle region
                 1.1149e+03 -5.0925e+03 0                                    % index 1 5 6 is a triangle region
                   985.6216 -5.0289e+03 0
                   865.9283 -5.1091e+03 0]';                                 % the order here is important and the order is counter-clock-wise
sample_points = Rx*sample_points;
sample_points(3,:) = atom_sim_pos(3,1);

%% set parameter
center_group = 2;                                                          % set as 1 or 2, needed be explained later
atomDist = 10;  

%% find the triangle region from the graphene system
unit3=[];
for tempi = 1:size(atom_sim_pos,2)
    if atom_group(tempi) == center_group
        % move the sample to the local region
        shift = atom_sim_pos(:,tempi) - sample_points(:,1);
        sample_points = sample_points + shift;

        % set the array to store the atom index match the sample points
        temp_array = zeros(1,6);
        for tempj = 1:size(sample_points,2)
            temp_dis = sqrt(sum((sample_points(:,tempj) - atom_sim_pos).^2,1));
            temp_ind = temp_dis<atomDist;
            temp_ind = find(temp_ind==1);
            if size(temp_ind,2) == 1
                temp_array(tempj) = temp_ind;
            end
            clear temp_ind temp_dis
        end
        clear tempj

        % group the triangle region from the graphene system
        if temp_array(2) ~= 0 && temp_array(3) ~=0
            % the index of the three vertices
            ind0 = [tempi temp_array(2) temp_array(3)];
            ind0 = sort(ind0);

            % calculate the average position center of the triangle (centroid)
            cen=mean(atom_sim_pos(:,ind0),2)';      
            unit3=[unit3;ind0,cen]; 
            clear cen ind0
        end
        if temp_array(3) ~= 0 && temp_array(4) ~=0
            % the index of the three vertices
            ind0 = [tempi temp_array(3) temp_array(4)];
            ind0 = sort(ind0);

            % calculate the average position center of the triangle (centroid)
            cen=mean(atom_sim_pos(:,ind0),2)';      
            unit3=[unit3;ind0,cen]; 
            clear cen ind0
        end
        if temp_array(4) ~= 0 && temp_array(5) ~=0
            % the index of the three vertices
            ind0 = [tempi temp_array(4) temp_array(5)];
            ind0 = sort(ind0);

            % calculate the average position center of the triangle (centroid)
            cen=mean(atom_sim_pos(:,ind0),2)';      
            unit3=[unit3;ind0,cen];
            clear cen ind0
        end
        if temp_array(5) ~= 0 && temp_array(6) ~=0
            % the index of the three vertices
            ind0 = [tempi temp_array(5) temp_array(6)];
            ind0 = sort(ind0);

            % calculate the average position center of the triangle (centroid)
            cen=mean(atom_sim_pos(:,ind0),2)';      
            unit3=[unit3;ind0,cen]; 
            clear cen ind0
        end
        clear shift temp_array
    end
end

% the code below delete the repeated row data 
unit3=unique(unit3,'rows');
clear sample_points tempi

%% find the neighbor triangular region
neighbor_ind = zeros(size(unit3,1),3);
for tempi = 1:size(neighbor_ind,1)
    ind_1 = unit3(tempi,1); 
    ind_2 = unit3(tempi,2);
    ind_3 = unit3(tempi,3);
    match_ind_1 = (unit3(:,1)==ind_1) + (unit3(:,2)==ind_1) + (unit3(:,3)==ind_1);
    match_ind_2 = (unit3(:,1)==ind_2) + (unit3(:,2)==ind_2) + (unit3(:,3)==ind_2);
    match_ind_3 = (unit3(:,1)==ind_3) + (unit3(:,2)==ind_3) + (unit3(:,3)==ind_3);
    match_ind = match_ind_1 + match_ind_2 + match_ind_3;
    match_ind = find(match_ind==2);
    for tempj = 1:size(match_ind,1)
        neighbor_ind(tempi,tempj) = match_ind(tempj);
    end
    clear ind_1 ind_2 ind_3
    clear match_ind match_ind_1 match_ind_2 match_ind_3
end

%% calculate the Meron Number
func_meron = @(a) solid_angle(a(1,:),a(2,:),a(3,:)) ./ (4*pi);                   % the function to calculate Meron Numb
Meron_number = zeros(size(unit3,1),4);
for tempi = 1:size(unit3,1)
    % obtain the position of the atom
    x=atom_sim_pos(1,unit3(tempi,1:3))';                                            % the x-axis positions of the three vertices                                         
    y=atom_sim_pos(2,unit3(tempi,1:3))';                                            % the y-axis positions of the three vertices
    z=atom_sim_pos(3,unit3(tempi,1:3))';                                            % the z-axis positions of the three vertices
    temp_ind = (boundary(x,y))';                                                 % obtain the boundary of the triangle and the direction of the boundary

    % obtain the displacement of current atom
    n1 = atom_displace(:,unit3(tempi,temp_ind(1)))';
    n2 = atom_displace(:,unit3(tempi,temp_ind(2)))';
    n3 = atom_displace(:,unit3(tempi,temp_ind(3)))';
    Meron_number(tempi,1:3) = unit3(tempi,4:end);
    Meron_number(tempi,4)   = func_meron([n1;n2;n3]);
    clear temp_tan temp_solid_angle 
    clear n1 n2 n3
end

X_Coordinates   = double(Meron_number(:,1));
Y_Coordinates   = double(Meron_number(:,2));
Z_Coordinates   = double(Meron_number(:,4));

x     = linspace(min(X_Coordinates),max(X_Coordinates),500);
dx    = abs(x(2)-x(1));
yspan = abs(max(Y_Coordinates) - min(Y_Coordinates));
y     = linspace(min(Y_Coordinates),max(Y_Coordinates),round(yspan/dx));
clear yspan

%% the first possible region
% Center of the circle
center = [1803.21525734349,3433.80149279143];
%% scan all possibility 
a_array = 750:50:850;
b_array = 1450:50:1550;
phi_array = 15:0.1:16;
% Number of points
numPoints = 1000;
% Angle from 0 to 2*pi
theta = linspace(0, 2*pi, numPoints);
% error value
error_matrix = zeros(size(a_array,2),size(b_array,2),size(phi_array,2));
error = Inf;

counta = 0;
for tempa = a_array
    tic
    counta = counta + 1;
    countb = 0;
    for tempb = b_array
        countb = countb + 1;
        countphi = 0;
        for tempphi = phi_array
            
            countphi = countphi + 1;
            % Parametric equations for the ellipse before rotation
            x_ellipse = tempa * cos(theta);
            y_ellipse = tempb * sin(theta);
            % Rotation matrix
            R = [cosd(tempphi) -sind(tempphi); sind(tempphi) cosd(tempphi)];
            % Parametric equations for the ellipse before rotation
            x_ellipse = tempa * cos(theta);
            y_ellipse = tempb * sin(theta);
            % Apply rotation and translate to center
            temp_points_1 = R * [x_ellipse; y_ellipse];
            temp_points_1(1,:) =   temp_points_1(1,:) + center(1);
            temp_points_1(2,:) = - temp_points_1(2,:) - center(2);

            isInside = inpolygon(atom_sim_pos(1,:), atom_sim_pos(2,:), temp_points_1(1,:), temp_points_1(2,:));
            temp_ind = find(isInside==1);
            temp_region_inf = zeros(size(unit3,1),1); 
            for tempi = temp_ind
                 temp_region_inf =  temp_region_inf + sum(unit3(:,1:3) == tempi,2); % store this vortex inside Meron information as a neighbor  
            end
            temp_Meron_number = sum(Meron_number(temp_region_inf==3,4));
            temp_error = abs(temp_Meron_number - 0.5);
            
            if temp_error < error
                error = temp_error;
                a = tempa;
                b = tempb;
                phi = tempphi;
                points_1 = temp_points_1;
            end

            error_matrix(counta,countb,countphi) = temp_error;
            clear x_ellipse y_ellipse R isInside temp_ind temp_region_inf tempi temp_Meron_number
        end
    end
    toc
end
clear error

isInside = inpolygon(atom_sim_pos(1,:), atom_sim_pos(2,:), points_1(1,:), points_1(2,:));
temp_ind = find(isInside==1);
temp_region_inf = zeros(size(unit3,1),1); 
for tempi = temp_ind
     temp_region_inf =  temp_region_inf + sum(unit3(:,1:3) == tempi,2);    % store this vortex inside Meron information as a neighbor  
end
temp_Meron_number = sum(Meron_number(temp_region_inf==3,4));
temp_ind = find(temp_region_inf==3);
temp_ind_stack = [];
for tempi = temp_ind
    temp_ind_stack = [temp_ind_stack neighbor_ind(tempi,:)];
end
temp_ind_stack(temp_ind_stack==0) =     [];
temp_ind_stack = unique(temp_ind_stack,'stable');
region_inf{1}.Meron_number = temp_Meron_number;  
region_inf{1}.index        = temp_ind;     
region_inf{1}.ind_stack    = reshape(temp_ind_stack,max(size(temp_ind_stack)),1); 
clear temp_ind temp_Meron_number temp_ind_stack isInside points_1

%% match the potential Meron or Anti-Meron to 0.5 or -0.5 Skyrmion or Anti-Skyrmion to 1 or -1
Meron = {};
anti_Meron = {};
Skyrmion = {};
anti_Skyrmion = {};
for Index = 1:length(region_inf)
    temp_Meron_num = region_inf{Index}.Meron_number;
    temp_ind       = region_inf{Index}.index';
    ind_stack      = region_inf{Index}.ind_stack';

    %% find the Meron Number Region that is equal to -0.5 Meron equal to 1 Skyrmion
    if temp_Meron_num < 0
        if abs(temp_Meron_num + 0.5) < 0.1 || temp_Meron_num >= -0.4
            ref_density = -0.5;
        elseif abs(temp_Meron_num + 1) < 0.1 || temp_Meron_num <= -0.6
            ref_density = -1;
        else
            ref_density = -0.5;
        end
        temp_neighbor_ind = neighbor_ind;
        temp_neighber_inf = zeros(size(neighbor_ind,1),2);                 % the first colume represent whether this index is inside the Meron the second colume represent how many neighbor of this vortex inside the Meron
        error = abs(temp_Meron_num-ref_density);                           % the error between reference density with current Meron number
        
        for tempi = temp_ind
            % store this vortex inside Meron information           
            temp_neighber_inf(tempi,1) = 1;        
            % store this vortex inside Meron information as a neighbor  
            temp_neighber_inf(:,2) =  temp_neighber_inf(:,2) + sum(neighbor_ind == tempi,2);
            temp_neighbor_ind(temp_neighbor_ind == tempi) = 0;             % to avoid repeat calculate              
            ind_stack(ind_stack==tempi) = [];                              % pop out from stack  
        end
        % loop for searching the Meron region
        while true          
            % end loop condition 
            if isempty(ind_stack)
                break
            elseif error < 1
                break
            end

            % find the region decrease the error most
            temp_error = error;
            stack_index = 1;
            for temp_stack_index = 1:size(ind_stack,2)
                if abs(temp_Meron_num+Meron_number(ind_stack(temp_stack_index),4)-ref_density) <= temp_error
                    temp_error = abs(temp_Meron_num+Meron_number(ind_stack(temp_stack_index),4)-ref_density);
                    stack_index = temp_stack_index;
                end
            end
            clear temp_stack_index
            if temp_error == error
                break
            end
            
            if temp_error <= error
                temp_ind(end+1) = ind_stack(stack_index);
                temp_neighber_inf(ind_stack(stack_index),1) = 1;                                                  % store this vortex inside Meron information                                                      
                temp_neighber_inf(:,2) =  temp_neighber_inf(:,2) + sum(neighbor_ind == ind_stack(stack_index),2); % store this vortex inside Meron information as a neighbor      
                ind_stack = [ind_stack temp_neighbor_ind(ind_stack(stack_index),:)];
                ind_stack = unique(ind_stack, 'stable');                                                          % to avoid repeating calculate
                temp_neighbor_ind(temp_neighbor_ind == ind_stack(stack_index)) = 0;                               % to avoid repeat calculate             
    
                temp_Meron_num = temp_Meron_num + Meron_number(ind_stack(stack_index),4);
                error = abs(temp_Meron_num-ref_density);
                ind_stack(stack_index) = [];                                                % pop out from stack         
            else
                ind_stack(stack_index) = [];                                                % pop out from stack         
            end
    
            temp_edge_region = find(temp_neighber_inf(:,1)==1 & temp_neighber_inf(:,2)<3)'; % store the edge region index
            temp_edge_points = zeros(3*length(temp_edge_region),1);                         % store the edge points
            count = 1;
            for tempj = temp_edge_region
                temp_edge_points(3*(count-1)+1:3*count) = unit3(tempj,1:3); 
                count = count + 1;
            end
            clear count
            temp_edge_points = unique(temp_edge_points);
            ind_stack(ind_stack==0) = []; 
    
            for tempj = ind_stack
                isInside = inpolygon(atom_sim_pos(1,unit3(tempj,1:3)),atom_sim_pos(2,unit3(tempj,1:3)),atom_sim_pos(1,temp_edge_points),atom_sim_pos(2,temp_edge_points));
                if sum(isInside) == 3
    
                    temp_ind(end+1) = tempj;
                    temp_neighber_inf(tempj,1) = 1;                                                  % store this vortex inside Meron information                                                      
                    temp_neighber_inf(:,2) =  temp_neighber_inf(:,2) + sum(neighbor_ind == tempj,2); % store this vortex inside Meron information as a neighbor
                    temp_neighbor_ind(temp_neighbor_ind == tempj) = 0;                               % to avoid repeat calculate
                    ind_stack = [ind_stack temp_neighbor_ind(tempj,:)];
                    ind_stack = unique(ind_stack, 'stable');               % to avoid repeating calculate
                    temp_neighbor_ind(temp_neighbor_ind == tempj) = 0;     % to avoid repeat calculate
    
                    temp_Meron_num = temp_Meron_num + Meron_number(tempj,4);
                    error = abs(temp_Meron_num-ref_density);
                    ind_stack(ind_stack==tempj) = [];                      % pop out from stack                       
                end
                clear isInside
            end
            ind_stack(ind_stack==0) = []; 
            clear temp_edge_region temp_edge_points temp_error stack_index
            
        end
        if error < 1
            temp_atom_ind = []; 
    
            figure(334+Index);hold on
            for tempj = temp_ind
                temp_atom_ind = [temp_atom_ind, unit3(tempj,1:3)];
            end
            temp_atom_ind = unique(temp_atom_ind);
    
            for tempj = temp_atom_ind
                quiver3(atom_sim_pos(1,tempj),atom_sim_pos(2,tempj),0, atom_displace(1,tempj)*100, atom_displace(2,tempj)*100, atom_displace(3,tempj)*100, 0, 'MaxHeadSize', 1.5, 'Color', 'r', 'LineWidth', 1.5);hold on;
            end
            title('Meron Number: ',num2str(temp_Meron_num))
            axis equal
            hold off;
            
            if abs(temp_Meron_num+0.5) == error
                Meron{end+1}.Meron_number = temp_Meron_num; 
                Meron{end}.region_index = temp_ind;  
                Meron{end}.atom_index   =temp_atom_ind; 

                % save the atom position and displacement
                atomup_sim = atom_sim_pos(:,temp_atom_ind);
                displacementUP = atom_displace(:,temp_atom_ind);
                region_index = temp_ind;
                atom_index = temp_atom_ind;
                save([pwd,'/output_data/BLG_upper_layer_meron_2.mat'],'atomup_sim','displacementUP','region_index','atom_index'); 
                clear atomup_sim displacementUP region_index atom_index

            elseif abs(temp_Meron_num+1) == error
                Skyrmion{end+1}.Meron_number = temp_Meron_num; 
                Skyrmion{end}.region_index = temp_ind;  
                Skyrmion{end}.atom_index   =temp_atom_ind;
                   
                % save the atom position and displacement
                atomup_sim = atom_sim_pos(:,temp_atom_ind);
                displacementUP = atom_displace(:,temp_atom_ind);
                region_index = temp_ind;
                atom_index = temp_atom_ind;
                save([pwd,'/output_data/BLG_upper_layer_skyrmion_2.mat'],'atomup_sim','displacementUP','region_index','atom_index');
                clear atomup_sim displacementUP region_index atom_index

            end

        end
        clear temp_neighbor_ind ind_stack 
    end
    
    %% find the anti-Meron Number Region that is equal to 0.5 Meron or equal to 1 anti-Skyrmion
    if temp_Meron_num > 0
        if abs(temp_Meron_num - 0.5) < 0.1 || temp_Meron_num <= 0.4
            ref_density = 0.5;
        elseif abs(temp_Meron_num - 1) < 0.1 || temp_Meron_num >= 0.6
            ref_density = 1;
        else
            ref_density = 0.5;
        end

        temp_neighbor_ind = neighbor_ind;
        temp_neighber_inf = zeros(size(neighbor_ind,1),2);                 % the first colume represent whether this index is inside the Meron the second colume represent how many neighbor of this vortex inside the Meron
        error = abs(temp_Meron_num-ref_density);                           % the error between reference density with current Meron number
        for tempi = temp_ind
            % store this vortex inside Meron information           
            temp_neighber_inf(tempi,1) = 1;        
            % store this vortex inside Meron information as a neighbor  
            temp_neighber_inf(:,2) =  temp_neighber_inf(:,2) + sum(neighbor_ind == tempi,2);
            temp_neighbor_ind(temp_neighbor_ind == tempi) = 0;             % to avoid repeat calculate              
            ind_stack(ind_stack==tempi) = [];                              % pop out from stack  
        end
        % loop for searching the Meron region
        while true            
            % end loop condition 
            if isempty(ind_stack)
                break
            elseif error < 1
                break
            end
            % find the region decrease the error most
            temp_error = error;
            stack_index = 1;
            for temp_stack_index = 1:size(ind_stack,2)
                if abs(temp_Meron_num+Meron_number(ind_stack(temp_stack_index),4)-ref_density) <= temp_error
                    temp_error = abs(temp_Meron_num+Meron_number(ind_stack(temp_stack_index),4)-ref_density);
                    stack_index = temp_stack_index;
                end
            end
            clear temp_stack_index
            if temp_error == error
                break
            end

            if temp_error <= error
                temp_ind(end+1) = ind_stack(stack_index);
                temp_neighber_inf(ind_stack(stack_index),1) = 1;                                                  % store this vortex inside Meron information                                                      
                temp_neighber_inf(:,2) =  temp_neighber_inf(:,2) + sum(neighbor_ind == ind_stack(stack_index),2); % store this vortex inside Meron information as a neighbor      
                ind_stack = [ind_stack temp_neighbor_ind(ind_stack(stack_index),:)];
                ind_stack = unique(ind_stack, 'stable');                                                          % to avoid repeating calculate
                temp_neighbor_ind(temp_neighbor_ind == ind_stack(stack_index)) = 0;                               % to avoid repeat calculate             
    
                temp_Meron_num = temp_Meron_num + Meron_number(ind_stack(stack_index),4);
                error = abs(temp_Meron_num-ref_density);
                ind_stack(stack_index) = [];                                                % pop out from stack         
            else
                ind_stack(stack_index) = [];                                                % pop out from stack         
            end
    
            temp_edge_region = find(temp_neighber_inf(:,1)==1 & temp_neighber_inf(:,2)<3)'; % store the edge region index
            temp_edge_points = zeros(3*length(temp_edge_region),1);                         % store the edge points
            count = 1;
            for tempj = temp_edge_region
                temp_edge_points(3*(count-1)+1:3*count) = unit3(tempj,1:3); 
                count = count + 1;
            end
            clear count
            temp_edge_points = unique(temp_edge_points);
            ind_stack(ind_stack==0) = []; 
    
            for tempj = ind_stack
                isInside = inpolygon(atom_sim_pos(1,unit3(tempj,1:3)),atom_sim_pos(2,unit3(tempj,1:3)),atom_sim_pos(1,temp_edge_points),atom_sim_pos(2,temp_edge_points));
                if sum(isInside) == 3
    
                    temp_ind(end+1) = tempj;
                    temp_neighber_inf(tempj,1) = 1;                                                  % store this vortex inside Meron information                                                      
                    temp_neighber_inf(:,2) =  temp_neighber_inf(:,2) + sum(neighbor_ind == tempj,2); % store this vortex inside Meron information as a neighbor
                    temp_neighbor_ind(temp_neighbor_ind == tempj) = 0;                               % to avoid repeat calculate
                    ind_stack = [ind_stack temp_neighbor_ind(tempj,:)];
                    ind_stack = unique(ind_stack, 'stable');                                         % to avoid repeating calculate
                    temp_neighbor_ind(temp_neighbor_ind == tempj) = 0;                               % to avoid repeat calculate
    
                    temp_Meron_num = temp_Meron_num + Meron_number(tempj,4);
                    error = abs(temp_Meron_num-ref_density);
                    ind_stack(ind_stack==tempj) = [];                                                % pop out from stack                       
                end 
                clear isInside
            end
            ind_stack(ind_stack==0) = []; 
            clear temp_edge_region temp_edge_points temp_error stack_index
        end
        if error < 1
            temp_atom_ind = []; 
    
            figure(334+Index);hold on
            for tempj = temp_ind
                temp_atom_ind = [temp_atom_ind, unit3(tempj,1:3)];
            end
            temp_atom_ind = unique(temp_atom_ind);
    
            for tempj = temp_atom_ind
                quiver3(atom_sim_pos(1,tempj),atom_sim_pos(2,tempj),0, atom_displace(1,tempj)*100, atom_displace(2,tempj)*100, atom_displace(3,tempj)*100, 0, 'MaxHeadSize', 1.5, 'Color', 'r', 'LineWidth', 1.5);hold on;
            end
            title('anti-Meron Number: ',num2str(temp_Meron_num))
            axis equal
            hold off;

            if abs(temp_Meron_num-0.5) == error
                anti_Meron{end+1}.Meron_number = temp_Meron_num; 
                anti_Meron{end}.region_index = temp_ind; 
                anti_Meron{end}.atom_index   = temp_atom_ind; 

                % save the atom position and displacement
                atomup_sim = atom_sim_pos(:,temp_atom_ind);
                displacementUP = atom_displace(:,temp_atom_ind);
                region_index = temp_ind;
                atom_index = temp_atom_ind;
                save([pwd,'/output_data/BLG_upper_layer_anti_meron_2.mat'],'atomup_sim','displacementUP','region_index','atom_index'); 
                clear atomup_sim displacementUP region_index atom_index

            elseif abs(temp_Meron_num-1) == error
                anti_Skyrmion{end+1}.Meron_number = temp_Meron_num; 
                anti_Skyrmion{end}.region_index = temp_ind; 
                anti_Skyrmion{end}.atom_index   = temp_atom_ind; 

                % save the atom position and displacement
                atomup_sim = atom_sim_pos(:,temp_atom_ind);
                displacementUP = atom_displace(:,temp_atom_ind);
                region_index = temp_ind;
                atom_index = temp_atom_ind;
                save([pwd,'/output_data/BLG_upper_layer_anti_skyrmion_2.mat'],'atomup_sim','displacementUP','region_index','atom_index'); 
                clear atomup_sim displacementUP region_index atom_index
            end

        end
        clear temp_neighbor_ind ind_stack
    end
    clear temp_Meron_num temp_neighber_inf temp_ind ind_stack
end




