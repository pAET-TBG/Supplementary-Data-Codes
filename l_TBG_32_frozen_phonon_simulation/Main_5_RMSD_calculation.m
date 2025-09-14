clear;
clc;

%% load the atom data
addpath([pwd, '/input_data/'])
addpath([pwd, '/output_data/'])
addpath([pwd, '/src/'])
atompos_1 = importdata('atom_tracing_model_refinement_total.mat');
atompos_2 = importdata('atom_tracing_model_original.mat');
atompos_1 = atompos_1./0.99495;
Res = 0.19/2;                                                              % resolution
atompos_2(1,:) = atompos_2(1,:) - mean(atompos_2(1,:)) + mean(atompos_1(1,:));
atompos_2(2,:) = atompos_2(2,:) - mean(atompos_2(2,:)) + mean(atompos_1(2,:));
atompos_2(3,:) = atompos_2(3,:) - mean(atompos_2(3,:)) + mean(atompos_1(3,:));
  
%% subpixel alignment
atomDist = 0.6;
Mag_shift=1;
shift=[0 0 0]';                                                            % the shift used to do the subpixel shift to align the atom position
while Mag_shift>1e-5                                                       % if the total difference is larger than 1e-6 (angstrom) then continue to do the alignment
    atompos_2 =  atompos_2 - shift;                                        % substruct the difference (i.e. do the shift to align the atom)
    difarr = [];                                                           % array for store the difference between two tracing result
    disarr = [];                                                           % array for store the distance between two tracing result
    count_arr1 = [];                                                       % store the index for the first atom result, which is match with the second atom tracing result
    count_arr2 = [];                                                       % store the index for the second atom result, which is match with the first atom tracing result
    for i=1:size(atompos_1,2)
        dif=(atompos_2-atompos_1(:,i));                                    % calculate all difference from the first result with i'th position in the second result (angstrom)
        dis=sqrt(sum(dif.^2,1));                                           % calculate the distance (angstrom)
        [dis,ind]=min(dis);                                                % obatin the minimum distance and corresponding index
        if dis <= atomDist                                                 % if the minimum distance is smaller than the threshold then store the information
            difarr=[difarr dif(:,ind)];
            disarr = [disarr dis];
            count_arr1 = [count_arr1,ind];
            count_arr2 = [count_arr2,i];
        end
    end
    shift=mean(difarr,2);                                                  % calculate the mean values of the difference in x-y-z axis, which is also the shift value for the next iteration
    Mag_shift=sum(abs(mean(difarr,2)));                                    % calculate the total difference
end

[dif12,disarr12]=calculate_model_difference(atompos_1,atompos_2,atomDist);
sig12=std(dif12,0,2);
meanD12 = mean(dif12,2);

figure(111); clf;
set(gcf,'position',[0,100,1000,600]);

subplot(131); histogram(abs(dif12(1,:)).*100,20);hold on;
title(['y-axis: RMSD: ',num2str(sig12(1).*100),'pm'],[' Mean: ',num2str(meanD12(1).*100), 'pm']);
xlabel('The position difference in x-axis (pm)')
ylabel('The number of atom')
set(gca,'fontsize',15,'linewidth',1);

subplot(132);histogram(abs(dif12(2,:)).*100,20);hold on;
title(['x-axis: RMSD: ',num2str(sig12(2).*100),'pm'],[' Mean: ',num2str(meanD12(2).*100), 'pm']);
xlabel('The position difference in y-axis (pm)')
ylabel('The number of atom')
set(gca,'fontsize',15,'linewidth',1);

subplot(133); histogram(abs(dif12(3,:)).*100,20);hold on;
title(['z-axis: RMSD: ',num2str(sig12(3).*100),'pm'],[' Mean: ',num2str(meanD12(3).*100), 'pm']);
xlabel('The position difference in z-axis (pm)')
ylabel('The number of atom')
set(gca,'fontsize',15,'linewidth',1);

%% calculate the total RMSD
rmsd_value12 = sqrt(sum(sig12.^2));
fprintf('rmsd = %.2f angstrom\n',rmsd_value12);
fprintf('rmsd = %.2f pm\n',rmsd_value12*100);

figure(12);
% Adjust its size
set(figure(12), 'Position', [1, 1, 500*1.8, 450*1.8]);hold on;
histogram(disarr12.*100,40,'FaceColor',[0 0.4470 0.7410],'EdgeColor', 'k','LineWidth',5/1.5,'BinWidth', 1.155/3*3.0);hold on;
xlim([0  30]);ylim([0 900]);                                               % Set the x-axis limits from 0 to 50
xticks(0:15:30);yticks(0:225:900);
% Increase the font size of the tick labels
ax = gca;                                                                  % Get the current axes
ax.Layer = 'top';                                                          % Draw the axes lines and labels over the histograms
ax.FontName = 'Arial';                                                     % Change font type of axis labels and ticks
ax.FontSize = 51/1.5;                                                      % Set the font size. Adjust this value as needed
% Make the axes box visible if it's not already
ax.Box = 'on';
% Set the line width of the box (frame line of the plot)
ax.LineWidth = 5/1.5;                                                      % Adjust this value to make the frame line thicker
hold off

%% function used to calculate the difference between two model
function [difarr, disarr] = calculate_model_difference(atompos_1,atompos_2,atomDist)
Mag_shift=1;
shift=[0 0 0]';                                                            % the shift used to do the subpixel shift to align the atom position
while Mag_shift>1e-5                                                       % if the total difference is larger than 1e-6 (angstrom) then continue to do the alignment
    atompos_2 =  atompos_2 - shift;                                        % substruct the difference (i.e. do the shift to align the atom)
    difarr = [];                                                           % array for store the difference between two tracing result
    disarr = [];                                                           % array for store the distance between two tracing result
    count_arr1 = [];                                                       % store the index for the first atom result, which is match with the second atom tracing result
    count_arr2 = [];                                                       % store the index for the second atom result, which is match with the first atom tracing result
    for i=1:size(atompos_1,2)
        dif=(atompos_2-atompos_1(:,i));                                    % calculate all difference from the first result with i'th position in the second result (angstrom)
        dis=sqrt(sum(dif.^2,1));                                           % calculate the distance (angstrom)
        [dis,ind]=min(dis);                                                % obatin the minimum distance and corresponding index
        if dis <= atomDist                                                 % if the minimum distance is smaller than the threshold then store the information
            difarr=[difarr dif(:,ind)];
            disarr = [disarr dis];
            count_arr1 = [count_arr1,ind];
            count_arr2 = [count_arr2,i];
        end
    end
    shift=mean(difarr,2);                                                  % calculate the mean values of the difference in x-y-z axis, which is also the shift value for the next iteration
    Mag_shift=sum(abs(mean(difarr,2)));                                    % calculate the total difference
end
end


