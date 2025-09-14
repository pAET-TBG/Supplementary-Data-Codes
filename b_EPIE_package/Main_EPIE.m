clear;
clc

addpath([pwd, '/src/']);
addpath([pwd, '/input_data/']);
addpath([pwd, '/output_data/']);
load('Angles_scan_id.mat')

for scan_id = 1:13
    load([pwd, '/input_data/diffraction_pattern_',num2str(scan_id),'.mat']);
    dp=m;
    scan_col=size(dp,3);
    scan_row=size(dp,4);
    dp_slices=sum(squeeze(sum(dp,3)),3);
    dp_back=dp_slices/(size(dp,3)*size(dp,4));
    for ii=1:size(dp,3)
        for jj=1:size(dp,4)
            dp_temp=dp(:,:,ii,jj);
            dp_temp(dp_temp<0)=0;
            dp(:,:,ii,jj) = dp_temp;
        end
    end
    dp_slices=sum(squeeze(sum(dp,3)),3);
    dp_sum= dp_slices;
    dp_slices(dp_slices<=1e9)=0;
    dp_slices(dp_slices>1e9)=1e9;
    [centers,radii] = imfindcircles(dp_slices,[10 70],'ObjectPolarity','bright', 'Sensitivity',0.95,'EdgeThreshold',0.1);
      
    dx_preset=0.21;
    N_set=size(dp,1)*4;
    N_set_center=round(N_set/2);
    dp_raw=zeros(N_set,N_set, scan_col, scan_row);
    dp_raw(N_set_center-round(centers(1,2))+1+1:(N_set_center+48-round(centers(1,2)))+1,N_set_center-round(centers(1,1))+1:(N_set_center+48-round(centers(1,1))),:,:)=dp;
    clear dp
    dp=single(dp_raw);
    N = size(dp,1);
    dp=reshape(dp, N, N, []);
    dp(dp<0)=0; 
    for ii=1:size(dp,3)
        for jj=1:size(dp,4)
            dp_temp=fftshift(dp(:,:,ii,jj));
            dp(:,:,ii,jj)=sqrt(dp_temp);
        end
    end
    dp_LQ=fftshift(gainDisc(N, (radii-4)*2));
    for ll = 1:size(dp,3)
        dp_temp=dp(:,:,ll);
        dp(:,:,ll) = dp_temp;
    end
    N = size(dp,1);
    
    %% parameters
    alpha_max =32;                                                         % aperture size(mrad), the convergence angle 11.5 mrad is the semi-angle
    voltage =70;                                                           % beam voltage (keV)
    scanStepSize_x = 0.652;                                                % scan step size in horizontal direction (A)
    scanStepSize_y = 0.652;                                                % scan step size in vertical direction (A)
    rotAngle = 90;                                                         % rotation angle between CBED and actual scan position (degree)
    lambda = 12.398/sqrt((2*511.0+voltage).*voltage);                      % angstrom
    dx =0.19/2;
        
    %% reconstruction
    fprintf('ePIE reconstruction...\n');
    Niter = 100;                                                           % total number of iteration
    [px, py]=calculateScanPositions(scan_row,scan_col , scanStepSize_x, scanStepSize_y, rotAngle);
    
    % reshape diffraction pattern and scan position
    px = px(:);
    py = py(:);
    dp_avg = sum(sum(dp,3),4)/size(dp,3)/size(dp,4);

    % normalize initial probe
    df=40;
    probe0=generateProbeFunction(dx,N,0,0,df,0,1,voltage,alpha_max,0);
    probe{1} =probe0;
    probe{2} = fftshift(imshift_fft(fftshift(probe0), 0, 0.5, true));
    probe{3} = fftshift(imshift_fft(fftshift(probe0), 0.5, 0, true));
    probe{4} = fftshift(imshift_fft(fftshift(probe0), 0, 0, true));
    probe_orth = ortho_modes_eig(probe);
    pr(:,:,1) =  probe_orth{1};
    pr(:,:,2) =  probe_orth{2};
    pr(:,:,3) =  probe_orth{3};
    pr(:,:,4) =  probe_orth{4};
    pr=pr*sum(sum(dp_avg))/sum(sum(sqrt(abs(fftn(pr(:,:,1))).^2+ abs(fftn(pr(:,:,2))).^2+ abs(fftn(pr(:,:,3))).^2 ...
        + abs(fftn(pr(:,:,4))).^2)));
    M=4;
    % initializae object function
    Ny_max = max(abs(round(min(py)/dx)-floor(N/2)), abs(round(max(py)/dx)+ceil(N/2)))*2+1;
    Nx_max = max(abs(round(min(px)/dx)-floor(N/2)), abs(round(max(px)/dx)+ceil(N/2)))*2+1;
    
    N_obj = max(Ny_max,Nx_max)+10;
    ind__obj_center = floor(N_obj/2)+1;
    obj_ini = ones(N_obj, N_obj);
    for nn = 1:1
        obj{nn}= obj_ini;
    end
    % calculate indicies for all scans
    N_scan = round(length(px)/1);
    % position = pi(integer) + pf(fraction)
    py_i = round(py/dx);
    py_i_p = py/dx;
    py_f = py - py_i*dx;
    px_i = round(px/dx);
    px_i_p = px/dx;
    px_f = px - px_i*dx;
    
    ind_x_lb = px_i - floor(N/2) + ind__obj_center;
    ind_x_ub = px_i + ceil(N/2) -1 + ind__obj_center;
    ind_y_lb = py_i - floor(N/2) + ind__obj_center;
    ind_y_ub = py_i + ceil(N/2) -1 + ind__obj_center;
    
    %update step size
    alpha = 0.3;
    beta = 0.3;
    alpha_g = 0.3;
    beta_g = 0.3;
    obj_modul=ones(N_obj);
    ob_shift = single(zeros(N,N,M));
    N_obj_l=round(N_obj*1.2);
    V=0-N_obj_l/2:0+N_obj_l/2-1;
    [X,Y]=meshgrid(V);
    Z=zeros(N_obj_l);
    
    phi   = -angles_pre(scan_id,1)/180*pi;                                 % z-axis
    theta = -angles_pre(scan_id,2)/180*pi;                                 % x-axis
    psi   = angles_pre(scan_id,3)/180*pi;                                  % y-axis

    [Xnew,Ynew,Znew]=AxisAngleRotate2(X,Y,Z,[1,1,1],[psi,theta,phi]);
    Xnew= round(Xnew);
    Ynew= round(Ynew);
    Xnew_Ynew=[Xnew(:),Ynew(:)];
    transform_posi=zeros(1,size(dp,3));
    for aper = 1:size(dp,3)
        ind_X=find(Xnew==px_i(aper));
        ind_Y=find(Ynew==py_i(aper));
        ind_XY=intersect(ind_X,ind_Y);
        if length(ind_XY)==0
            transform_posi(aper)=transform_posi(aper-1);
        else
            transform_posi(aper)=Znew(ind_XY(1))*dx;
        end
    end

    cm=0.01;um=1e-6;                                     
    lambda_t=lambda*0.0001*um;
    MM=N;                              
    hx  = N*dx*0.0001*um; hy  = N*dx*0.0001*um;                               
    dhx = hx / MM;        dhy = hy / N;
    k   = 2*pi / lambda_t;
    du  = 1 ./ (MM*dhx);
    dv  = 1 ./ (N*dhy);
    u   = ones(N,1)*[0:MM/2-1 -MM/2+1:0]*du;                               % Note order of points for FFT
    v   = [0:N/2-1 -N/2+1:0]'*ones(1,MM)*dv;
       
    %% position correction
    usfac       =1000;
    int_dp=squeeze(sum(sum(dp)));
    int_dp_beta_possion=(int_dp-min(int_dp(:)))/(max(int_dp(:))-min(int_dp(:)));
    int_dp_beta_possion(int_dp_beta_possion<0.1)=0.01;
    int_dp_beta_possion(int_dp_beta_possion>=0.1)=0.3;
    probe_region=zeros(N);
    for ii=1:N
        for jj=1:N
            if sqrt((ii-N/2)^2+(jj-N/2)^2)<25 
                probe_region(ii,jj)=1;
            end
        end
    end
    offset_corr_x=-(px_i_p-px_i);
    offset_corr_y=-(py_i_p-py_i);
    
    [nr,nc]=size(dp(:,:,1));
    Nr = ifftshift(-fix(nr/2):ceil(nr/2)-1);
    Nc = ifftshift(-fix(nc/2):ceil(nc/2)-1);
    [Nc,Nr] = meshgrid(Nc,Nr);
    gpu=1;
    if gpu==1
        obj = cellfun(@gpuArray, obj, 'UniformOutput', false);
        pr = gpuArray(pr);
        obj_modul=gpuArray(obj_modul);
        ob_shift = gpuArray(ob_shift);
        probe_region=  gpuArray(probe_region);
    end
    
    posi_rand_x=round(rand(N_scan,1)*4-2);
    posi_rand_y=round(rand(N_scan,1)*4-2);

    ind_y_lb_rand=ind_y_lb;
    ind_y_ub_rand=ind_y_ub;
    ind_x_lb_rand=ind_x_lb;
    ind_x_ub_rand=ind_x_ub;
    for iter=1:100
        fprintf(strcat('Iteration:',num2str(iter),'/',num2str(Niter),'\n'));
        updateOrder = randperm(N_scan);
        
        error_psi=0;
        for ii=1:N_scan
            ind_dp = updateOrder(ii);
            if  ismember(ind_dp,[])
            else  
                ob0_shift_main = obj{1}(ind_y_lb_rand(ind_dp):ind_y_ub_rand(ind_dp),ind_x_lb_rand(ind_dp):ind_x_ub_rand(ind_dp));
                dp0=dp(:,:,ind_dp);
                dp1 = zeros(N,N);
                for mm = 1:M
                    if mm<=4
                        pr0 = pr(:,:,mm);
                        H=exp(1i*2*pi*(0.0001*transform_posi(ind_dp)*um)*sqrt((1/lambda_t).^2-u.^2-v.^2));        
                        pr0=ifft2(fft2(pr0).*H); 
                        pr_temp(:,:,mm)=pr0;
                        psi0 = pr0.*ob0_shift_main;
                        Psi0 = fftn(psi0);
                        ob_shift(:,:,mm) = ob0_shift_main;
                        psi_all(:,:,mm) = psi0;
                        Psi_all(:,:,mm) = Psi0;
                        dp1 = dp1 + abs(Psi0).^2;
                    end   
                end
                dp1 = sqrt(dp1);
                dp_up = dp0./(dp1+eps);
                if ind_dp>=1&&ind_dp<=225
                    error_psi=error_psi+sum(sum(abs(dp0-dp1)));
                end
                err_array(iter,1)=gather(error_psi);
                
                pr0_acc_main = zeros(N,N);
                pr0_diff_acc_main = zeros(N,N);
                for mm = 1:M
                    ob0_shift = ob_shift(:,:,mm);
                    pr0 = pr_temp(:,:,mm);
                    psi0 = psi_all(:,:,mm);
                    Psi0 = Psi_all(:,:,mm);
                    Psi1 = dp_up .* Psi0;
                    psi1 = ifftn(Psi1);
                    if mm<=4
                        if iter<=0
                            pr1 = pr0;
                        else
                            pr1 = pr0 + alpha *conj(ob0_shift)/(max(abs(ob0_shift(:)))^2).*(psi1 - psi0);
                        end
                        pr0_acc_main=  pr0_acc_main+abs(pr0).^2;
                        pr0_diff_acc_main = pr0_diff_acc_main+ beta *conj(pr0).*(psi1 - psi0);
                    end
                    H=exp(1i*2*pi*(0.0001*-transform_posi(ind_dp)*um)*sqrt((1/lambda_t).^2-u.^2-v.^2));        
                    pr1=ifft2(fft2(pr1).*H); 
                    pr(:,:,mm) = pr1;
                end
                obj_illumed_pre = probe_region.*ob0_shift_main;
                ob1_shift_main = ob0_shift_main + pr0_diff_acc_main/max(pr0_acc_main(:));
                obj_illumed_now = probe_region.*ob1_shift_main;
                obj{1}(ind_y_lb_rand(ind_dp):ind_y_ub_rand(ind_dp),ind_x_lb_rand(ind_dp):ind_x_ub_rand(ind_dp)) = ob1_shift_main;
            end
            pr=pr*sum(sum(dp_avg))/sum(sum(sqrt(abs(fftn(pr(:,:,1))).^2+ abs(fftn(pr(:,:,2))).^2+ abs(fftn(pr(:,:,3))).^2 + abs(fftn(pr(:,:,4))).^2)));
        end
        
        if mod(iter,100)==0 
            obj_main=gather(obj{1});
            save([pwd,'/output_data/EPIE_reconstruction_',num2str(scan_id),'.mat'],'obj_main');

            pr=gather(pr);
            save([pwd,'/output_data/EPIE_probe_',num2str(scan_id),'.mat'],'pr','err_array');
        end
    end
end