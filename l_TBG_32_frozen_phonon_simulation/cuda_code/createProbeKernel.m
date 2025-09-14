function probe_Kernel = createProbeKernel(V, pixelsize, nr, alpha, df_series)%#codegen
%% parameters
h=6.62*1.0e-34;% Planck's constant
e=1.6*1.0e-19;% charge of electron
m=9.1*1.0e-31; % mass of electron, kg
%V=300*10^3;% 300 keV
c=3*1.0e8; % speed of light
lamuda=h/sqrt(2*m*e*V)*1/sqrt(1+e*V/2/m/c^2); % electron wave length in m
lamuda=lamuda*10^10; % electron wave length in Angst.
%nr=200; % number of pixel of the probe
%pixelsize=0.332; % real space pixelsize in Angst.
%alpha=17.1; % convergence semi-angle
%df_series=100;% Defocus
%% construct aperture
kmax=1/pixelsize*lamuda/2*1e3; % mrad
dkx=(kmax-(-kmax))/nr;

kx=(-kmax):dkx:kmax;
[kX, kY]=meshgrid(kx, kx);

%kX_rad=kX*0.001;% mrad
%kY_rad=kY*0.001;% mrad
%kR_rad = sqrt(kX.^2+kY.^2);% mrad
kR_rad = kX.^2+kY.^2;% mrad

ape = kR_rad;
ape(ape<=alpha^2)=1;
ape(ape>alpha^2)=0;
%figure(2); img(ape)

%% construct probe


%w = kX_rad + kY_rad*1i;
%w_conj = kX_rad - kY_rad*1i;
%base_func0 = -w.*w_conj/2;
%base_func  = -(kX_rad.^2 + kY_rad.^2 )/2;
%base_func  = -(kX.^2 + kY.^2 )*(1e-6/2);
%norm(base_func0-base_func,'fro')/norm(base_func0,'fro')

%abe_func = df_series*base_func*2*pi/lamuda;
abe_func = kR_rad* (-df_series * 1e-6*pi/lamuda);

probe_func = exp(-1i*abe_func);
probe_func = probe_func.*ape;
probe_func = ifft2(ifftshift(probe_func));
probe_func = (fftshift(probe_func));
%probe_func_conj = conj(probe_func);
%probe_Kernel0 = probe_func.*probe_func_conj;
probe_Kernel  = abs(probe_func).^2;
%norm(probe_Kernel0-probe_Kernel,'fro')/norm(probe_Kernel0,'fro')
%figure; img(probe_Kernel);

%% Normalize
sum_Probe = sum(probe_Kernel,'all');
probe_Kernel = probe_Kernel./sum_Probe;
end



