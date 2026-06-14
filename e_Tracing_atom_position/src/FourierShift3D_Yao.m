function shift_vol = FourierShift3D_Yao(vol,dx,dy,dz)
[dimx_pj, dimy_pj, dimz_pj] = size(vol);
[Y_pj, X_pj, Z_pj] = meshgrid(-(dimx_pj-1)/2:(dimx_pj-1)/2,...
    -(dimy_pj-1)/2:(dimy_pj-1)/2,-(dimz_pj-1)/2:(dimz_pj-1)/2); 

shift_mat = exp((2j*pi/dimy_pj*dy)*Y_pj +...
    (2j*pi/dimx_pj*dx)* X_pj +(2j*pi/dimz_pj*dz)* Z_pj);
% shift_mat = exp((-2j*pi/dimy_pj*dy)*Y_pj +...
%     (-2j*pi/dimx_pj*dx)* X_pj +(-2j*pi/dimz_pj*dz)* Z_pj);
vol_fft = my_fft(vol);
shift_vol = my_ifft(vol_fft.*shift_mat);
end