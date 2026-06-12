clear;
clc;

%% load data
addpath([pwd,'\src\'])
addpath([pwd,'\input_data\'])
addpath([pwd,'\output_data\'])

%% generate simulation atom model
atom_sim = generate_simulation_model();
atom_sim = atom_sim.*0.9145;
atom_sim(1,:) = atom_sim(1,:) + 0.2;
atom_sim(2,:) = atom_sim(2,:) - 50.4;
%% check the refinement coordination
atompos_1 = importdata('atom_tracing_model_refinement_total_GPU.mat');
atompos_1 = atompos_1./0.99495;
points    = [-31 -49
              31 -49  
              31  49
             -31  49]'; 
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
atomDist = 1;
label_atom  = zeros(1,size(atom_sim,2));
for tempi = 1:size(atompos_1,2)
    dif=(atom_sim-atompos_1(:,tempi));         % calculate all difference from the first result with i'th position in the second result (angstrom)
    dis=sqrt(sum(dif.^2,1));                  % calculate the distance (angstrom)
    [dis,ind]=min(dis);                       % obatin the minimum distance and corresponding index
    if dis < atomDist                        % if the minimum distance is smaller than the threshold then store the information
        label_atom(:,ind) = 1;
    end
end
atom_sim = atom_sim(:,label_atom == 1);

% match the atom position for the first round
[atompos_1,atom_sim] = atom_align_exp_sim_model_1(atompos_1,atom_sim);
atom_exp = atompos_1;
[atom_exp,atom_sim] = atom_align_exp_sim_model_1(atom_exp,atom_sim);
[atom_exp,atom_sim] = atom_align_exp_sim_model_2(atom_exp,atom_sim);
[atom_exp,atom_sim] = atom_align_exp_sim_model_3(atom_exp,atom_sim);
[atom_exp,atom_sim] = atom_align_exp_sim_model_4(atom_exp,atom_sim);
[atom_exp,atom_sim] = atom_align_exp_sim_model_5(atom_exp,atom_sim);

%% match each experiment atom with simulation atom
atom = zeros(size(atom_exp));
for tempi = 1:size(atom_exp,2)
    dif=(atom_exp(:,tempi)-atom_sim);                                          % calculate all difference from the first result with i'th position in the second result (angstrom)
    dis=sqrt(sum(dif.^2,1));                                           % calculate the distance (angstrom)
    [dis,ind]=min(dis);  
    if dis <= 0.5  % obatin the minimum distance and corresponding index
        atom(:,tempi) = atom_sim(:,ind);
    end
end
atom_sim = atom;
atomup = atom_exp(:,atom_exp(3,:)<0).*100;
atomdown = atom_exp(:,atom_exp(3,:)>0).*100;
atomup_sim = atom_sim(:,atom_sim(3,:)<0).*100;
atomdown_sim = atom_sim(:,atom_sim(3,:)>0).*100;

save([pwd,'/output_data/atom_tracing_model_refinement_total_matched.mat'],'atomup','atomdown','atom_exp','atomup_sim','atomdown_sim','atom_sim')

