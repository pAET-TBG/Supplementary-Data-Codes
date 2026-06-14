function [temp_model, temp_atomtype,score] = initial_classification1_1types_fixnonatom(...
    RecVol_padded, curr_model, classify_info,bg)
% Yao Yang

% 2024.04.10
% fix non-atom shape equal to zero
% 2019.05.16
% from class7 to classify two type

% two step: first classify non-atoms from the originally method
% (use three type, which already modified)
% then select certain number of atoms around 3 peaks, calculate the
% averaged atom size then move them freely

% Modified with Minh's implement:
% 1. reduced several sub-functions
% 2. add L1 and L2 norm options
% Notice: curr_model should be (3 by N)

lnorm = 2;

if isfield(classify_info,'phase_shift_flag')       phase_shift_flag = classify_info.phase_shift_flag; %#ok<SEPEX>
else phase_shift_flag = 1;          %#ok<SEPEX>
end
if isfield(classify_info,'halfSize')       halfSize = classify_info.halfSize; %#ok<SEPEX>
else halfSize = 1;          %#ok<SEPEX>
end
if isfield(classify_info,'plothalfSize')   plothalfSize = classify_info.plothalfSize; %#ok<SEPEX>
else plothalfSize = 4;      %#ok<SEPEX>
end
if isfield(classify_info,'SPHyn')          SPHyn = classify_info.SPHyn; %#ok<SEPEX>
else SPHyn = 1;             %#ok<SEPEX>
end
if isfield(classify_info,'separate_part')  separate_part = classify_info.separate_part; %#ok<SEPEX>
else separate_part = 100;   %#ok<SEPEX>
end
if isfield(classify_info,'L_forAver')      L_forAver = classify_info.L_forAver; %#ok<SEPEX>
else L_forAver = 100;       %#ok<SEPEX>
end
if isfield(classify_info,'boundary')       boundary = classify_info.boundary;
else boundary = 5;
end
% curr_model: 3 by n
if phase_shift_flag == 1
    compute_average_func = @compute_average_atom_from_vol_fouriershift;
else
    compute_average_func = @compute_average_atom_from_vol;
end
%%
dim = size(RecVol_padded,1);
b1 = find(curr_model(1,:)<boundary | curr_model(1,:)>dim-boundary+1);
b2 = find(curr_model(2,:)<boundary | curr_model(2,:)>dim-boundary+1);
b3 = find(curr_model(3,:)<boundary | curr_model(3,:)>dim-boundary+1);
bT = union(union(b1,b2),b3);
curr_model(:,bT) = [];
%% obtain global intensity histogram

[xX,yY,zZ] = ndgrid(-halfSize:halfSize,-halfSize:halfSize,-halfSize:halfSize);
SphereInd = find(xX.^2+yY.^2+zZ.^2 <=(halfSize+0.5)^2);
[xXp,yYp,zZp] = ndgrid(-plothalfSize:plothalfSize,-plothalfSize:plothalfSize,-plothalfSize:plothalfSize);
SphereInd_plot = find(xXp.^2+yYp.^2+zZp.^2 <=(plothalfSize+0.5)^2);

if SPHyn
  useInd = SphereInd;
  useInd_plot =SphereInd_plot;
else
  useInd = 1:length(xX);
  useInd_plot = 1:length(xXp);
end

%% compute initial average Fe atom, Pt atom and non atom

% re-produce array for integrated intensity
intensity_integ = zeros(1,size(curr_model,2));
intensity_integ_plot = zeros(1,size(curr_model,2));

% integrate intensity (for 3x3x3 and 5x5x5 voxels) for each traced peak
for j=1:size(curr_model,2)
    curr_pos  = round(curr_model(:,j));
    box_integ = RecVol_padded(curr_pos(1)-halfSize:curr_pos(1)+halfSize,...
        curr_pos(2)-halfSize:curr_pos(2)+halfSize,curr_pos(3)-halfSize:curr_pos(3)+halfSize);
    intensity_integ(j) = sum(box_integ(useInd));
end
for j=1:size(curr_model,2)
    curr_pos  = round(curr_model(:,j));
    box_integ = RecVol_padded(curr_pos(1)-plothalfSize:curr_pos(1)+plothalfSize,...
        curr_pos(2)-plothalfSize:curr_pos(2)+plothalfSize,curr_pos(3)-plothalfSize:curr_pos(3)+plothalfSize);
    intensity_integ_plot(j) = sum(box_integ(useInd_plot));
end

%%
[hist_integ, cen_integ] = hist(intensity_integ_plot,separate_part);
initcen_integ   = [cen_integ(round(separate_part/2))/2, cen_integ(round(separate_part/2))*3/2];
initpeak_integ  = [max(hist_integ) max(hist_integ)];
initwidth_integ = [cen_integ(round(separate_part/2))/4, cen_integ(round(separate_part/2))/4]*3;
i_guess = [0 initpeak_integ(1) initcen_integ(1) initwidth_integ(1) 0 initpeak_integ(2) initcen_integ(2) initwidth_integ(2)];
Xdata = cen_integ;
Ydata = hist_integ; 
% Ydata(1:20) = 0;
lb = [-Inf -Inf 0 -Inf -Inf -Inf -Inf -Inf];
ub = [ Inf  Inf 0  Inf  Inf  Inf  Inf  Inf];

[p, ~, fitresult] = My_two_gaussianfit(Xdata, Ydata, i_guess);

% figure(101)
% hist(intensity_integ_plot,separate_part);
% hold on
% plot(cen_integ, fitresult, 'r-', 'LineWidth',3);
% plot(cen_integ, (p(1)+p(5))/2+p(2)*exp(-((Xdata-p(3))/p(4)).^2),'g-','LineWidth',2);
% plot(cen_integ, (p(1)+p(5))/2+p(6)*exp(-((Xdata-p(7))/p(8)).^2),'g-','LineWidth',2);
% hold off
% xlabel('integrated intensity (a.u.)');
% ylabel('# atoms');
% title(sprintf('pos [%5.2f, %5.2f]\n width [%5.2f, %5.2f] height [%5.2f, %5.2f]',...
%     p(3), p(7), abs(p(4)), abs(p(8)), p(2), p(6)));
%%
% curr_model(:,intensity_integ_plot<50000) = [];
% intensity_integ(intensity_integ_plot<50000) = [];
% intensity_integ_plot(intensity_integ_plot<50000) = [];
atomtype = zeros(1,length(intensity_integ));

%%
avg_atom = [];

[sort_intensity_integ,sort_ind] = sort(intensity_integ_plot);

[~,minInd_1] = min(abs(sort_intensity_integ-p(3)));
lower_bound_1 = max(1,minInd_1-L_forAver);
if minInd_1+L_forAver <= size(sort_ind,2)
    def_1_ind = sort_ind(lower_bound_1:minInd_1+L_forAver);
else
    def_1_ind = sort_ind(lower_bound_1:size(sort_ind,2));
end
[avg_atom_temp]= compute_average_func(RecVol_padded,curr_model,def_1_ind,halfSize);
avg_atom = [avg_atom,avg_atom_temp(useInd)];

[~,minInd_2] = min(abs(sort_intensity_integ-p(7)));
lower_bound_2 = max(1,minInd_2-L_forAver);
upper_bound_2 = min(length(intensity_integ_plot),minInd_2+L_forAver);
def_2_ind = sort_ind(lower_bound_2:upper_bound_2);
[avg_atom_temp]= compute_average_func(RecVol_padded,curr_model,def_2_ind,halfSize);
avg_atom = [avg_atom,avg_atom_temp(useInd)];

atomtype(def_1_ind) = 1;
atomtype(def_2_ind) = 2;
%% generate points of intensities 
Num_atom = size(curr_model,2);
pt_size  = length(useInd);
points   = zeros(pt_size, Num_atom);
for k = 1:size(curr_model,2)
    curr_x = round(curr_model(1,k));
    curr_y = round(curr_model(2,k));
    curr_z = round(curr_model(3,k));
    
    if phase_shift_flag == 1
        dx = curr_model(1,k) - curr_x;
        dy = curr_model(2,k) - curr_y;
        dz = curr_model(3,k) - curr_z;
    end
    curr_vol = RecVol_padded(curr_x-halfSize-1:curr_x+halfSize+1,...
        curr_y-halfSize-1:curr_y+halfSize+1,...
        curr_z-halfSize-1:curr_z+halfSize+1);
    if phase_shift_flag == 1
        shiftModel = FourierShift3D_Yao(curr_vol,dx,dy,dz);
        shift_vol = shiftModel(2:end-1,2:end-1,2:end-1);
    else
        shift_vol = curr_vol(2:end-1,2:end-1,2:end-1);
    end
    
    points(:,k) = shift_vol(useInd);
end

%% atom classification iteration of Types
dist_mat = zeros(2,1);

%intensity_integ_plot = intensity_integ_plot.^2;
intensity_integ_plot = intensity_integ_plot-min(intensity_integ_plot);
avg_atom(:,2) = avg_atom(:,2)-min(points(:));
points = points-min(points(:));
%bg=0.015
for iter=1:101
    
    % fix non-atom shape equal to 0
   bg*mean(intensity_integ_plot);

    avg_atom(:,1) = bg*mean(intensity_integ_plot);

    [hist_inten_plot,cen_integ_total_plot]= hist(intensity_integ_plot,separate_part);
    y_up = max([round(max(hist_inten_plot)/10)*12, 5]);
    
    figure(103)
    clf
    subplot(3,1,1);
    hist(intensity_integ_plot,separate_part);
    %xlim([0 ceil(max(intensity_integ_plot)/5)*5]);
    xlim([0 (max(intensity_integ_plot))]);
    ylim([0 y_up]);
    xlabel('integrated intensity (a.u.)');
    ylabel('# atoms');
    title(sprintf('boxsize %d',halfSize*2+1));
    
    for i = 1:2
        intensity_integ_sub = intensity_integ_plot(atomtype==i);
        subplot(3,1,i+1)
        
        hist(intensity_integ_sub,cen_integ_total_plot);
        xlabel('integrated intensity (a.u.)');
        ylabel('# atoms');
        title(sprintf('%d Type %i atoms',sum(atomtype==i),i));
        ylim([0 y_up]);
    %    xlim([0 ceil(max(intensity_integ_plot)/5)*5]);
        xlim([0 (max(intensity_integ_plot))]);
    end
    drawnow;
%     pause()
    old_labels = atomtype;

    obj = 0;
    for n = 1:Num_atom
        for k=1:2
            dist_mat(k) = norm(points(:,n) - avg_atom(:,k), lnorm).^lnorm;            
        end
        [dist2, idx] = min(dist_mat);        
%        [dist2, idx] = min( sum( bsxfun(@minus, avg_atom, points(:,n) ).^2)   );
        atomtype(n) = idx;
        obj = obj + dist2;
    end
    for k=1:2
        avg_atom(:,k) = sum(points(:,atomtype==k),2)/sum(atomtype==k);
    end
    
    score=(obj/Num_atom)^(1/lnorm);

    fprintf('%02i. obj = %.3f\n',iter, (obj/Num_atom)^(1/lnorm) );
    
    % if there is no change in the atomic specise classification, break
    if ~any(old_labels-atomtype), break; end
    
end
temp_model = curr_model;
temp_atomtype = atomtype;
end