clear;
clc;

%% load data
addpath([pwd,'/src/'])
addpath([pwd,'/input_data/'])

%% generate simulation atom model
atom_sim = generate_simulation_model();

%% the parameter to judge the quality of each traced atom
Res = 0.2055/2;                    
% rotation matrix
theta0 = -55.2;
R0 = [cosd(theta0) -sind(theta0); sind(theta0) cosd(theta0)];
atom_sim(1:2,:) = R0*(atom_sim(1:2,:));
atom_sim = atom_sim./Res;
atom_sim(1,:) = atom_sim(1,:)+173;
atom_sim(2,:) = atom_sim(2,:)+562;
atom_sim(3,:) = atom_sim(3,:)+40;

%% check the refinement coordination
load('atom_tracing_model_refinement.mat')
Res = 0.19/2;                    
atom_exp(1,:) = atom(1,:)./Res+824;
atom_exp(2,:) = atom(2,:)./Res+722;
atom_exp(3,:) = atom(3,:)./Res+40;
Index = (atom_sim(1,:)>=min(atom_exp(1,:))-400) & (atom_sim(1,:)<=max(atom_exp(1,:))+400) & (atom_sim(2,:)>=min(atom_exp(2,:))-250) & (atom_sim(2,:)<=max(atom_exp(2,:))+250);
atom_sim = atom_sim(:,Index);
clear atom m n R0 R1 R2 theta0 Index

%% center the coordinates
atom_exp(1,:) = atom_exp(1,:) - mean(atom_sim(1,:)) + 6.0949;
atom_exp(2,:) = atom_exp(2,:) - mean(atom_sim(2,:)) - 10.8449;
atom_exp(3,:) = atom_exp(3,:) - mean(atom_sim(3,:)) + 0.0784;
atom_sim(1,:) = atom_sim(1,:) - mean(atom_sim(1,:)) + 6.0949;
atom_sim(2,:) = atom_sim(2,:) - mean(atom_sim(2,:)) - 10.8449;
atom_sim(3,:) = atom_sim(3,:) - mean(atom_sim(3,:)) + 0.0784;
atom_sim = atom_sim.*Res;
atom_exp = atom_exp.*Res;
clear Res

% match the atom position for the first round
[atom_exp,atom_sim] = atom_align_exp_sim_model_1(atom_exp,atom_sim);

%% check the twisting angle
theta0 = 1.9;
atomDist = 70;
searchRange = -0.1:0.001:0.1;
twist_angle = calculate_twist_angle(atom_exp, theta0, searchRange,atomDist, false);
fprintf('The twist angle: %.2f\n', twist_angle);
clear theta0 atomDist searchRange twist_angle

%% match each experiment atom with simulation atom
load('atom_tracing_model_refinement.mat');
atom     = atom.*100;
atom_exp = atom_exp.*100;

atom_exp(1,:) = atom_exp(1,:) - mean(atom(1,:));
atom_exp(2,:) = atom_exp(2,:) - mean(atom(2,:));
atom_exp(3,:) = atom_exp(3,:) - mean(atom(3,:));
atom(1,:) = atom(1,:) - mean(atom(1,:));
atom(2,:) = atom(2,:) - mean(atom(2,:));
atom(3,:) = atom(3,:) - mean(atom(3,:));

%% match the atom to find the label for Carbon and Silicon for simulation data
label_exp = zeros(size(label));
atomDist = 50;
Mag_shift=1;
shift=[0 0 0]';                                                            % the shift used to do the subpixel shift to align the atom position
while Mag_shift>1e-4                                                       % if the total difference is larger than 1e-4 (angstrom) then continue to do the alignment
    atom_exp =  atom_exp - shift;                                          % substruct the difference (i.e. do the shift to align the atom)
    difarr = [];                                                           % array for store the difference between two tracing result
    disarr = [];                                                           % array for store the distance between two tracing result
    for i=1:size(atom,2)
        dif=(atom_exp-atom(:,i));                                          % calculate all difference from the first result with i'th position in the second result (angstrom)
        dis=sqrt(sum(dif.^2,1));                                           % calculate the distance (angstrom)
        [dis,ind]=min(dis);                                                % obatin the minimum distance and corresponding index
        if dis <= atomDist                                                 % if the minimum distance is smaller than the threshold then store the information
            difarr=[difarr dif(:,ind)];
            disarr = [disarr dis];
            label_exp(i) = label(ind); 
        end
    end
    shift=mean(difarr,2);                                                  % calculate the mean values of the difference in x-y-z axis, which is also the shift value for the next iteration
    Mag_shift=sum(abs(mean(difarr,2)));                                    % calculate the total difference
end
clear Mag_shift shift dis disarr dif difarr i ind
clear atom label atomDist

atomup_exp    = atom_exp(:,atom_exp(3,:)<0);                               % store the inital position for distance calculate
atomup_sim    = importdata('atom_tracing_model_upper_sim.mat');            % load the corresponding simulation atom model
labelup       = label_exp(atom_exp(3,:)<0);                                % store the original label for each atom
atomupgroup = importdata('atom_tracing_model_upper_group.mat');
atomdown_exp  = atom_exp(:,atom_exp(3,:)>0);                               % store the inital position for distance calculate
atomdown_sim  = importdata('atom_tracing_model_lower_sim.mat');            % load the corresponding simulation atom model
labeldown     = label_exp(atom_exp(3,:)>0);                                % store the original label for each atom
atomdowngroup = importdata('atom_tracing_model_lower_group.mat');        
[atomup,atomup_sim,labelup,atomdown,atomdown_sim,labeldown,atomupgroup,atomdowngroup] =  crop_region(atomup_exp,atomdown_exp,atomup_sim,atomdown_sim,labelup,labeldown,atomupgroup,atomdowngroup);

%% match the crop region
[atomup,atomup_sim,atomdown,atomdown_sim] = atom_align_exp_sim_model_2(atomup,atomup_sim,atomdown,atomdown_sim);
save([pwd,'/output_data/atom_model_sim_exp_matched.mat'],'atomup','atomup_sim','labelup','atomupgroup','atomdown','atomdown_sim','labeldown','atomdowngroup')

