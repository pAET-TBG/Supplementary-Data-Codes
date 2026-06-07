clear;clc;
addpath([pwd,'\src\'])
addpath([pwd,'\input_data\'])

%% load the ptychography projections
Projs = zeros(1781,1781,13);
index = 1:13;
for tempi = index
    temp_projs = importdata(['EPIE_reconstruction_',num2str(tempi),'.mat']);  
    Projs(:,:,tempi) = rot90(angle(temp_projs) - min(min(angle(temp_projs))),-1);
end

projections(:,:,7) = Projs(:,:,1);   
projections(:,:,2) = Projs(:,:,2);  
projections(:,:,1) = Projs(:,:,3); 
projections(:,:,8) = Projs(:,:,4); 
projections(:,:,9) = Projs(:,:,5);   
projections(:,:,5) = Projs(:,:,6);   
projections(:,:,6) = Projs(:,:,7);   
projections(:,:,3) = Projs(:,:,8);   
projections(:,:,11) = Projs(:,:,9);  
projections(:,:,12) = Projs(:,:,10);
projections(:,:,10) = Projs(:,:,11); 
projections(:,:,4) = Projs(:,:,12);  
projections(:,:,13) = Projs(:,:,13); 
clear Projs
Projs = projections;

%% use 2D Gaussian to find the peak
[X,Y] = meshgrid(-5:6,-5:6); 
xdata.x = X;
xdata.y = Y;
fixed = [0 0 0 0 0 0 0];
lb = [-10 -10 -9 -9 1e1 1e1 -pi];
ub = [10 10 10 10 1e3 1e3 pi];
init_guess = [0 3 0 0 3e2 3e2 0];
index = [1 2 3 4 5 6 7 8 9 10 11 12 14];
inital_pos = [863 972; 931 991; 956 907; 965 1069;
              940 957; 1009 909; 997 991; 906 994;
              1026 1065; 786 967; 930 1123; 996 1024;
              873 864; 933 990; 1010 1039; 1051 765];
inital_pos = inital_pos(index,:);
fit_pos = size(inital_pos);

for iter = 1:5
    for tempi = 1:size(Projs,3)
        
        % plot the image to check whether the marker point is right
        figure(111);
        imagesc(Projs(:,:,tempi));
        hold on;
        
        % find the local region
        ydata = Projs(inital_pos(tempi,1)-5:inital_pos(tempi,1)+6, ...
                         inital_pos(tempi,2)-5:inital_pos(tempi,2)+6,tempi);
    
        % fit with 2D Gaussian
        [xfit, resnorm,residual] = fit_gauss2D_PD(init_guess, xdata, ydata, fixed, lb,ub);
        fit_pos(tempi,:) = inital_pos(tempi,:)+[xfit(3),xfit(4)];
        
        % scatter the fit point
        scatter(fit_pos(tempi,2),fit_pos(tempi,1),'r.');
        hold off
        drawnow;
        
        % catch whether the point is right
%         yfit = calc_gauss2D_PD(xfit,xdata);
%         figure();img(yfit,[],ydata,[])
    end
    inital_pos = round(fit_pos);
end
%% recenter the projections
center_projs = zeros(2500,2500,size(Projs,3));
for tempi = 1:size(Projs,3)
    col = round(fit_pos(tempi,1));
    row = round(fit_pos(tempi,2));
    % background substruction
    center_projs(1250-col+1:1250-col+size(Projs,1), ...
                 1250-row+1:1250-row+size(Projs,2),tempi) = Projs(:,:,tempi) - mean(mean(Projs(:,:,tempi)));
end
center_projs(:,:,1)  = imtranslate(center_projs(:,:,1),[190, 68]);
center_projs(:,:,2)  = imtranslate(center_projs(:,:,2),[188, 60]);
center_projs(:,:,3)  = imtranslate(center_projs(:,:,3),[193, 52]);
center_projs(:,:,4)  = imtranslate(center_projs(:,:,4),[192, 51]);
center_projs(:,:,5)  = imtranslate(center_projs(:,:,5),[188, 68]);
center_projs(:,:,6)  = imtranslate(center_projs(:,:,6),[181, 53]);
center_projs(:,:,7)  = imtranslate(center_projs(:,:,7),[188, 63]);
center_projs(:,:,8)  = imtranslate(center_projs(:,:,8),[185, 68]);
center_projs(:,:,9)  = imtranslate(center_projs(:,:,9),[181, 59]);
center_projs(:,:,10)  = imtranslate(center_projs(:,:,10),[176, 65]);
center_projs(:,:,11)  = imtranslate(center_projs(:,:,11),[182, 78]);
center_projs(:,:,12)  = imtranslate(center_projs(:,:,12),[194, 63]);
center_projs(:,:,13)  = imtranslate(center_projs(:,:,13),[190, 61]);
center_projs = center_projs - min(center_projs(:));
clear Projs;
Projs = center_projs(1250-999:1250+1000,1250-999:1250+1000,:);

%% normalization
for tempi = 1:size(Projs,3)
    Projs(:,:,tempi) = Projs(:,:,tempi) - min(min(Projs(1000-499:1000+500,1000-499:1000+500,tempi)));
    Projs(:,:,tempi) = Projs(:,:,tempi) ./max(max(Projs(1000-499:1000+500,1000-499:1000+500,tempi)));
    Projs(:,:,tempi) = Projs(:,:,tempi) - mean(mean(Projs(1000-499:1000+500,1000-499:1000+500,tempi)));
end
Projs = Projs - min(min(min(Projs(1000-499:1000+500,1000-499:1000+500,:))));
save([pwd,'\output_data\Projections_Rough_Alignment_Normalized.mat'],'Projs')



