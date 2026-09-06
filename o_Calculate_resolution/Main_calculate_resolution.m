clear;
clc;

%% load data
addpath([pwd,'/src/'])
addpath([pwd,'/input_data/'])
load('atom_tracing_model_refinement.mat')
load('tBLG_reconstruction_volume_without_support.mat')
rec = rec(1+200:end-200,1+200:end-200,:);
clear rec_1 mask_1 mask

atom_fit = importdata('atom_reference_volume.mat');
rec = rec(270:1364,65:1514,:);
atom = atom./(0.19/2);
atom(1,:) = atom(1,:) - mean(atom(1,:)) + size(rec,1)./2 + 24 - 17;
atom(2,:) = atom(2,:) - mean(atom(2,:)) + size(rec,2)./2 - 60 + 9;
atom(3,:) = atom(3,:) - mean(atom(3,:)) + size(rec,3)./2;

%% calculate the cross correlation
atom_mask = atom_fit;
atom_mask(atom_mask<0.0026) = nan;
atom_mask(atom_mask>0) = 1;
corr_vol = zeros(size(rec));
rec = My_paddzero(rec,[size(rec,1)+10,size(rec,2)+10,size(rec,3)+34]);
Ny = size(corr_vol,1);
Nx = size(corr_vol,2);
Nz = size(corr_vol,3);

parfor tempy = 1:Ny
    tic
    posy = 5+tempy;
    for tempx = 1:Nx
        posx = 5+tempx;
        for tempz = 1:Nz
            posz = 17+tempz;
            atom_temp = rec(posy-5:posy+5,posx-5:posx+5,posz-17:posz+17);
            valid_idx = ~isnan(atom_temp.*atom_mask) & ~isnan(atom_fit.*atom_mask) & isfinite(atom_temp.*atom_mask) & isfinite(atom_fit.*atom_mask);
            corr_vol(tempy,tempx,tempz) = corr(atom_temp(valid_idx), atom_fit(valid_idx));
        end
    end
    toc
end
clear Ny Nx Nz posy posx posxz atom_temp

%% Located the good atom
max_pos = [];
max_val = [];
count = 0;
for tempi = 1:size(atom,2)
    if label(tempi) == 1
        tempx = round(atom(2,tempi));
        tempy = round(atom(1,tempi));
        tempz = round(atom(3,tempi));
        temp_vol = corr_vol(tempy-3:tempy+3,tempx-3:tempx+3,tempz-3:tempz+3);
        if max(temp_vol(:)) >= 0.69
            [max_cor, linear_idx] = max(temp_vol(:));              
            [y, x, z] = ind2sub(size(temp_vol), linear_idx);
            if tempz - 4 + z >= 18 && tempz - 4 + z <= 63
                count = count+1;
                max_pos(:,count) = [tempy - 4 + y;tempx - 4 + x;tempz - 4 + z];
                max_val(count)   = max_cor;
            end
            clear z y x max_cor linear_idx
        end
        clear tempx tempy tempz temp_vol
    end
end

%% re-load the reconstruction
load('tBLG_reconstruction_volume_without_support.mat')
rec = rec(1+200:end-200,1+200:end-200,:);
rec = rec(270:1364,65:1514,:);
rec_pad = My_paddzero(rec,size(rec)+[0 0 8]);
rec_large = My_paddzero(rec,size(rec)+[10 10 38]);

%% Average the good atom
atom_tot = zeros(size(atom_fit)+[0 0 0]);
atom_tot_large  = zeros(23,23,67);
for tempi = 1:size(max_pos,2)
    tic
    comp_corr = 0;
    dx = 0;
    dy = 0;
    for tempy = -2:2
        for tempx = -2:2
            atom_temp = rec(max_pos(1,tempi)-5+tempy:max_pos(1,tempi)+5+tempy,max_pos(2,tempi)-5+tempx:max_pos(2,tempi)+5+tempx,max_pos(3,tempi)-17:max_pos(3,tempi)+17);
            valid_idx = ~isnan(atom_temp.*atom_mask) & ~isnan(atom_fit.*atom_mask) & isfinite(atom_temp.*atom_mask) & isfinite(atom_fit.*atom_mask);
            temp_corr = corr(atom_temp(valid_idx), atom_fit(valid_idx));
            if temp_corr > comp_corr
                comp_corr = temp_corr;
                dx = tempx;
                dy = tempy;
            end
            clear temp_corr atom_temp valid_idx
        end
    end
    clear tempx tempy comp_corr 

    temp_vol = rec_pad(max_pos(1,tempi)-5+dy:max_pos(1,tempi)+5+dy,max_pos(2,tempi)-5+dx:max_pos(2,tempi)+5+dx,max_pos(3,tempi)-21+4:max_pos(3,tempi)+21+4);
    comp_corr = 0;
    dz = 0;
    for tempz = -3:3
        atom_temp = temp_vol(:,:,22-17+tempz:22+17+tempz);
        valid_idx = ~isnan(atom_temp.*atom_mask) & ~isnan(atom_fit.*atom_mask) & isfinite(atom_temp.*atom_mask) & isfinite(atom_fit.*atom_mask);
        temp_corr = corr(atom_temp(valid_idx), atom_fit(valid_idx));
        if temp_corr > comp_corr
            comp_corr = temp_corr;
            dz = tempz;
        end
        clear temp_corr atom_temp valid_idx
    end
    clear tempz  
    clear atom_temp comp_corr 
    
    % subpixel alignment
    temp_vol = rec_pad(max_pos(1,tempi)-6+dy:max_pos(1,tempi)+6+dy,max_pos(2,tempi)-6+dx:max_pos(2,tempi)+6+dx,max_pos(3,tempi)-18+dz+4:max_pos(3,tempi)+18+dz+4);
    temp_vol_large = rec_large(max_pos(1,tempi)+5-6-5+dy:max_pos(1,tempi)+5+6+5+dy,max_pos(2,tempi)+5-6-5+dx:max_pos(2,tempi)+5+6+5+dx,max_pos(3,tempi)+15-18-15+dz+4:max_pos(3,tempi)+15+18+15+dz+4);
    corr_arr = zeros([21,21,21]);
    
    parfor tempy = 1:21
        for tempx = 1:21
            for tempz = 1:21
                atom_temp = abs(My_FourierShift_3D(temp_vol,(tempy-11)/10,(tempx-11)/10,(tempz-11)/10));
                atom_temp = atom_temp(7-5:7+5,7-5:7+5,19-17:19+17);
                valid_idx = ~isnan(atom_temp.*atom_mask) & ~isnan(atom_fit.*atom_mask) & isfinite(atom_temp.*atom_mask) & isfinite(atom_fit.*atom_mask);
                corr_arr(tempy, tempx, tempz) = corr(atom_temp(valid_idx), atom_fit(valid_idx));
            end
        end 
    end
    clear atom_temp valid_idx
    clear tempx tempy tempz  
    [max_cor, linear_idx] = max(abs(corr_arr(:)));     

    [dy, dx, dz] = ind2sub(size(corr_arr), linear_idx);
    clear linear_idx
    atom_temp = abs(My_FourierShift_3D(temp_vol,(dy-11)/10,(dx-11)/10,(dz-11)/10));
    atom_tot = atom_tot+abs(atom_temp(7-5:7+5,7-5:7+5,19-17:19+17));
    atom_tot_large = atom_tot_large + abs(My_FourierShift_3D(temp_vol_large,(dy-11)/10,(dx-11)/10,(dz-11)/10));
    toc
end
atom_tot = atom_tot./count;
atom_tot_large = atom_tot_large./size(max_pos,2);

%% Calculate the resolution
% normalized the averaged atom
atom_tot_large = atom_tot_large(:,:,2:end-1);
atom_tot_large = atom_tot_large - min(atom_tot_large(:));
atom_tot_large = atom_tot_large ./max(atom_tot_large(:));
% for y-direction
ydata = squeeze(squeeze(atom_tot_large(:,12,22)))';
xdata = 1:length(ydata);
init_guess = [ min(ydata), max(ydata)-min(ydata), xdata(ydata == max(ydata)), 30];
for iter = 1:5
    [x_fit_1D_y, ~, ~] = fit_gauss1D_PD(init_guess, xdata, ydata);
    init_guess = x_fit_1D_y;
end
% for x-direction
ydata = squeeze(atom_tot_large(12,:,22));
xdata = 1:length(ydata);
init_guess = [ min(ydata), max(ydata)-min(ydata), xdata(ydata == max(ydata)), 30];
for iter = 1:5
    [x_fit_1D_x, ~, ~] = fit_gauss1D_PD(init_guess, xdata, ydata);
    init_guess = x_fit_1D_x;
end
% for z-direction
ydata = squeeze(atom_tot_large(12,12,:))';
xdata = 1:length(ydata);
init_guess = [ min(ydata), max(ydata)-min(ydata), xdata(ydata == max(ydata)), 30];
for iter = 1:5
    [x_fit_1D_z, ~, ~] = fit_gauss1D_PD(init_guess, xdata, ydata);
    init_guess = x_fit_1D_z;
end

fprintf(['y-direction sigma ',num2str(sqrt(x_fit_1D_y(4)/2)*(0.19/2)),' A \n'])
fprintf(['x-direction sigma ',num2str(sqrt(x_fit_1D_x(4)/2)*(0.19/2)),' A \n'])
fprintf(['z-direction sigma ',num2str(sqrt(x_fit_1D_z(4)/2)*(0.19/2)),' A \n'])

fprintf(['y-direction resolution ',num2str(sqrt(x_fit_1D_y(4)/2)*(0.19/2)*2*sqrt(2*log(2))),' A \n'])
fprintf(['x-direction resolution ',num2str(sqrt(x_fit_1D_x(4)/2)*(0.19/2)*2*sqrt(2*log(2))),' A \n'])
fprintf(['z-direction resolution ',num2str(sqrt(x_fit_1D_z(4)/2)*(0.19/2)*2*sqrt(2*log(2))),' A \n'])

