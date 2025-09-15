function vol2 = My_FourierShift_3D(vol, dy, dx, dz)
% My_FourierShift_3D: Shifts a 3D volume in (dy, dx, dz) using Fourier shift theorem
%
% Inputs:
%   vol - 3D volume to shift
%   dy  - shift along Y (1st dim)
%   dx  - shift along X (2nd dim)
%   dz  - shift along Z (3rd dim)
%
% Output:
%   vol2 - shifted volume (same size)

% Get dimensions
[ny, nx, nz] = size(vol);

% Create centered frequency grids
[Y, X, Z] = ndgrid(-ceil((ny-1)/2):floor((ny-1)/2), ...
                   -ceil((nx-1)/2):floor((nx-1)/2), ...
                   -ceil((nz-1)/2):floor((nz-1)/2));

% Forward FFT
F = My_IFFTN(vol);  % User-defined inverse FFT

% Fourier phase shift factor
Pfactor = exp(2*pi*1i * (dx*X/nx + dy*Y/ny + dz*Z/nz));

% Apply shift and inverse FFT
vol2 = My_FFTN(F .* Pfactor);  % User-defined forward FFT
end