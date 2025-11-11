function [temp_model, temp_atomtype] = initial_class_kmean_sub_plot_ColorMap(rec, curr_model, classify_info,y_up)
% Yao Yang UCLA, 2020.08.05
% n-type classification by k-mean
% only classify among real-atoms

if isfield(classify_info,'lnorm')            lnorm = classify_info.lnorm; %#ok<SEPEX>
else lnorm = 2;             %#ok<SEPEX>
end
if isfield(classify_info,'Num_species')      Num_species = classify_info.Num_species; %#ok<SEPEX>
else Num_species = 3;       %#ok<SEPEX>
end
if isfield(classify_info,'halfSize')         halfSize = classify_info.halfSize; %#ok<SEPEX>
else halfSize = 1;          %#ok<SEPEX>
end
if isfield(classify_info,'plothalfSize')     plothalfSize = classify_info.plothalfSize; %#ok<SEPEX>
else plothalfSize = 4;      %#ok<SEPEX>
end
if isfield(classify_info,'separate_part')    separate_part = classify_info.separate_part; %#ok<SEPEX>
else separate_part = 70;      %#ok<SEPEX>
end
if isfield(classify_info,'O_Ratio')          O_Ratio = classify_info.O_Ratio; %#ok<SEPEX>
else O_Ratio = 1;           %#ok<SEPEX>
end
if isfield(classify_info,'SPHyn')            SPHyn = classify_info.SPHyn; %#ok<SEPEX>
else SPHyn = 1;             %#ok<SEPEX>
end
if isfield(classify_info,'PLOT_YN')          PLOT_YN = classify_info.PLOT_YN; %#ok<SEPEX>
else PLOT_YN = 0;           %#ok<SEPEX>
end
%% generate points of intensities 

[box_inten] = get_box_intensity(rec, curr_model, halfSize, O_Ratio, SPHyn, 'linear');
[box_inten_plot] = get_box_intensity(rec, curr_model, plothalfSize, O_Ratio, SPHyn, 'linear');

box_inten_plot=box_inten_plot/((2*halfSize)^3);
%% k-mean main function
if lnorm == 2
    idx = kmeans(sum(box_inten,1)',Num_species,'Start','uniform');
%     idx = kmeans(points',Num_species,'Distance','sqeuclidean');
elseif lnorm == 1
    idx = kmeans(sum(box_inten,1)',Num_species,'Distance','cityblock','Start','uniform');
end
%% Alignment clustered type into correct species order
mean_arr = zeros(1,Num_species);
for i = 1:Num_species
    mean_arr(i) = mean(sum(box_inten(:,idx==i),1));
end
[~,sortMean] = sort(mean_arr);
for i = 1:Num_species
    idx(idx==sortMean(i)) = i+1000;
end
idx = idx-1000;

%% final results
temp_model = curr_model;
temp_atomtype = idx';
%% plot histogram
Total_inten_plot=sum(box_inten_plot,1);
[hist_inten_plot,cen_integ_total_plot] = hist(Total_inten_plot,separate_part);
sprintf('Total No. of atoms: %d', length((idx)))


%y_up = round(max(hist_inten_plot)*1.1/3);

bin=cen_integ_total_plot(2)-cen_integ_total_plot(1);
Edges_integ_total_plot=cen_integ_total_plot-bin/2;

%colors = hsv(100);

colors = [220 111 0
    0 0 160]/255;

figure; clf;box on;
set(gcf,'Position',[150,150,350,220])
hold on
for i=1:Num_species
intensity_integ_sub = sum(box_inten_plot(:,(temp_atomtype==i)),1);
(sprintf('%d Type %i atoms',sum(idx==i),i))
color_id = randi([1, 100]);    % Random unique indices
selected_colors = colors(i, :);  % Extract the colors

FaceColor=selected_colors;
histogram(intensity_integ_sub, 'BinEdges', Edges_integ_total_plot, 'EdgeColor', 'black', ...
'FaceColor', FaceColor, 'FaceAlpha', 0.5,'linewidth',1);

lgd = legend('Non-atoms','C atoms');
lgd.Box = 'off'; 
end

set(gca,'Layer','top')
h = findobj(gca,'Type','patch');
set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1);
xlabel('Integrated intensity (a.u.)');
ylabel('Number of atoms');
%ylim([0 y_up]);
%ylim([0 3000]);
%ylim([0 3500]);
%ylim([0 4000]);
ylim([0 y_up]);
xlim([0 max(sum(box_inten_plot,1))*1.00]);

end