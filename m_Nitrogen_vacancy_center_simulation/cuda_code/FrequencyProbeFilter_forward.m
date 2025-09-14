function A_blurry = FrequencyProbeFilter_forward(A, kR_rad2, ape, df_series, lambda , hsize)
padding = 'circular';
sizeA = size(A);
num_slices = numel(df_series);
assert(num_slices==sizeA(3));

A_blurry = zeros(sizeA, 'single');

% disp(size(A))
% disp(hsize) 
% disp(padding)
%figure(1); img(A)

B = padImage(A, hsize, padding);
fftSize = size(B);
B = fft2(B);

for k = 1:num_slices
    %% probe function
    abe_func     = kR_rad2* (-df_series(k) * 1e-6*pi/lambda);

    probe_func   = exp(-1i*abe_func);
    probe_func   = probe_func.*ape;
    probe_func   = ifft2(ifftshift(probe_func));
    probe_func   = (fftshift(probe_func));
    probe_Kernel = abs(probe_func).^2;
    sum_Probe    = sum(probe_Kernel,'all');
    probe_Kernel = probe_Kernel'./sum_Probe;

    %% convolution

    B_i = ifft2( B(:,:,k) .* fft2(probe_Kernel, fftSize(1), fftSize(2)), 'symmetric' );

    A_blurry(:,:,k) = unpadImage(B_i, [sizeA(1),sizeA(2)]);

end
    %figure; img(A_blurry)
end



function padSize = computePadSize(sizeA, sizeH)

rankA = numel(sizeA);
rankH = numel(sizeH);

sizeH = [sizeH ones(1,rankA-rankH)];

padSize = floor(sizeH/2);

end



function [A, padSize] = padImage(A, hsize, padding)
padSize = computePadSize(size(A), hsize);
if ischar(padding)
    method = padding;
    padVal = [];
else
    method = 'constant';
    padVal = padding;
end
A = padarray_algo(A, padSize, method, padVal, 'both');
end



function A = unpadImage(A, outSize)
start = 1 + size(A) - outSize;
stop  = start + outSize - 1;
subCrop.type = '()';
subCrop.subs = {start(1):stop(1), start(2):stop(2)};
for dims = 3 : ndims(A)
    subCrop.subs{dims} = start(dims):stop(dims);
end
A = subsref(A, subCrop);

end