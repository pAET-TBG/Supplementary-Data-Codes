clear;
clc;

%% load data
addpath([pwd,'\src\'])
addpath([pwd,'\input_data\'])

%% generate simulation atom model
atom_sim = generate_simulation_model();
theta0 = -55.2;
R0 = [cosd(theta0) -sind(theta0); sind(theta0) cosd(theta0)];
atom_sim(1:2,:) = R0*(atom_sim(1:2,:));
atom_sim = atom_sim.*0.9245;
atom_sim(1,:) = atom_sim(1,:) - 0.1;
atom_sim(2,:) = atom_sim(2,:) - 39.9;
%% check the refinement coordination
atompos_1 = importdata('atom_tracing_model_refinement_total_GPU.mat');
atompos_1 = atompos_1.*1.0053;
points    = [-30 -42
             -21 -60
              14 -59
              35 -22
              32  26
              18  50
             -14  39
             -30   4]'; 
isInside_1  = inpolygon(atompos_1(1,:), atompos_1(2,:), points(1,:), points(2,:));
atompos_1   = atompos_1(:,isInside_1);
label_atom  = zeros(1,size(atompos_1,2));
atomDist = 1;
for tempi = 1:size(atompos_1,2)
    dif=(atompos_1-atompos_1(:,tempi));         % calculate all difference from the first result with i'th position in the second result (angstrom)
    dis=sqrt(sum(dif.^2,1));                  % calculate the distance (angstrom)
    [dis,ind]=min(dis);                       % obatin the minimum distance and corresponding index
    if dis < atomDist                        % if the minimum distance is smaller than the threshold then store the information
        label_atom(:,ind) = 1;
    end
end
atompos_1 = atompos_1(:,label_atom == 1);
clear label_atom
atomDist = 1.4;
label_atom  = zeros(1,size(atom_sim,2));
for tempi = 1:size(atompos_1,2)
    dif=(atom_sim-atompos_1(:,tempi));         % calculate all difference from the first result with i'th position in the second result (angstrom)
    dis=sqrt(sum(dif.^2,1));                  % calculate the distance (angstrom)
    label_atom(:,dis<=atomDist) = 1;
end
atom_sim = atom_sim(:,label_atom == 1);

%% group the atom
label_group = zeros(1,size(atompos_1,2));
label_group([1 4]) = 1;
inside_1 = [1 4];
inside_2 = [];
while sum(label_group==0) > 0
    for tempi = inside_1
        dif=(atompos_1-atompos_1(:,tempi));                    % calculate all difference from the first result with i'th position in the second result (angstrom)
        dis=sqrt(sum(dif.^2,1));    
        inside_2 = [inside_2 find(dis>0&dis<1.9&label_group==0)];
        label_group(dis>0&dis<1.9) = 2;
        inside_1(1) = [];
    end    
    for tempi = inside_2
        dif=(atompos_1-atompos_1(:,tempi));                    % calculate all difference from the first result with i'th position in the second result (angstrom)
        dis=sqrt(sum(dif.^2,1));    
        inside_1 = [inside_1 find(dis>0&dis<1.9&label_group==0)];
        label_group(dis>0&dis<1.9) = 1;
        inside_2(1) = [];
    end
end

label_group_sim = zeros(1,size(atom_sim,2));
label_group_sim([1 2654]) = 2;
inside_2 = [1 3401];
inside_1 = [];
while sum(label_group_sim==0) > 0    
    for tempi = inside_2
        dif=(atom_sim-atom_sim(:,tempi));                    % calculate all difference from the first result with i'th position in the second result (angstrom)
        dis=sqrt(sum(dif.^2,1));    
        inside_1 = [inside_1 find(dis>0&dis<1.9&label_group_sim==0)];
        label_group_sim(dis>0&dis<1.9) = 1;
        inside_2(1) = [];
    end
    for tempi = inside_1
        dif=(atom_sim-atom_sim(:,tempi));                    % calculate all difference from the first result with i'th position in the second result (angstrom)
        dis=sqrt(sum(dif.^2,1));    
        inside_2 = [inside_2 find(dis>0&dis<1.9&label_group_sim==0)];
        label_group_sim(dis>0&dis<1.9) = 2;
        inside_1(1) = [];
    end
end

atomDist = 0.6;
Mag_shift=1;
shift=[0 0 0]';   % the shift used to do the subpixel shift to align the atom position
while Mag_shift>1e-5 % if the total difference is larger than 1e-6 (angstrom) then continue to do the alignment
    atom_sim =  atom_sim - shift; % substruct the difference (i.e. do the shift to align the atom)
    difarr = [];    % array for store the difference between two tracing result
    disarr = [];    % array for store the distance between two tracing result
    count_arr1 = [];% store the index for the first atom result, which is match with the second atom tracing result
    count_arr2 = [];% store the index for the second atom result, which is match with the first atom tracing result
    for i=1:size(atompos_1,2)
        dif=(atom_sim-atompos_1(:,i));         % calculate all difference from the first result with i'th position in the second result (angstrom)
        dis=sqrt(sum(dif.^2,1));                  % calculate the distance (angstrom)
        [dis,ind]=min(dis);                       % obatin the minimum distance and corresponding index
        if dis <= atomDist                        % if the minimum distance is smaller than the threshold then store the information
            difarr=[difarr dif(:,ind)];
            disarr = [disarr dis];
            count_arr1 = [count_arr1,ind];
            count_arr2 = [count_arr2,i];
        end
    end
    shift=mean(difarr,2);                  % calculate the mean values of the difference in x-y-z axis, which is also the shift value for the next iteration
    Mag_shift=sum(abs(mean(difarr,2)));          % calculate the total difference
end
atomDist = 1.4;
label_atom  = zeros(1,size(atom_sim,2));
atom_sim_new = zeros(size(atompos_1));
for tempi = 1:size(atompos_1,2)
    if label_group(tempi) == 1
        index_sim = find(label_group_sim==1);
    elseif label_group(tempi) == 2
        index_sim = find(label_group_sim==2);
    end
    dif=(atom_sim(:,index_sim)-atompos_1(:,tempi));         % calculate all difference from the first result with i'th position in the second result (angstrom)
    dis=sqrt(sum(dif.^2,1));                  % calculate the distance (angstrom)
    [dis,ind]=min(dis);                       % obatin the minimum distance and corresponding index
    if dis < atomDist                        % if the minimum distance is smaller than the threshold then store the information
        atom_sim_new(:,tempi) = atom_sim(:,index_sim(ind));
    end
end
atom_sim = atom_sim_new;

% match the atom position for the first round
[atompos_1,atom_sim] = atom_align_exp_sim_model_1(atompos_1,atom_sim);
atom_exp = atompos_1;
[atom_exp,atom_sim] = atom_align_exp_sim_model_1(atom_exp,atom_sim);
[atom_exp,atom_sim] = atom_align_exp_sim_model_2(atom_exp,atom_sim);
[atom_exp,atom_sim] = atom_align_exp_sim_model_3(atom_exp,atom_sim);
[atom_exp,atom_sim] = atom_align_exp_sim_model_4(atom_exp,atom_sim);
[atom_exp,atom_sim] = atom_align_exp_sim_model_5(atom_exp,atom_sim);

atomup = atom_exp(:,atom_exp(3,:)<0).*100;
atomdown = atom_exp(:,atom_exp(3,:)>0).*100;
atomup_sim = atom_sim(:,atom_sim(3,:)<0).*100;
atomdown_sim = atom_sim(:,atom_sim(3,:)>0).*100;

save([pwd,'/output_data/atom_tracing_model_refinement_total_matched.mat'],'atomup','atomdown','atom_exp','atomup_sim','atomdown_sim','atom_sim')

