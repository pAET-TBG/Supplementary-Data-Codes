function TraceVol = Sub_Traced_Atom(volume,atompos,vol_size,iter)
if iter <= 6
    sigma_filter = 5;
elseif iter >= 16
    sigma_filter = 6;
else
    sigma_filter = 0.1*iter + 4.4;
end
atom_vol = zeros(2*vol_size+1,2*vol_size+1,2*vol_size+1);
filter = create3DGeneralizedGaussianFilter(vol_size,sigma_filter,2);
filter = filter./max(filter(:)).*1;
count = 0;
for tempi = 1:size(atompos,2)
    if ((round(atompos(2,tempi))-vol_size)>0) && ((round(atompos(1,tempi))-vol_size)>0) && ((round(atompos(3,tempi))-vol_size)>0) ...
       && (round(atompos(2,tempi))+vol_size<size(volume,2)) && (round(atompos(1,tempi))+vol_size<size(volume,1)) && (round(atompos(3,tempi))+vol_size<size(volume,3))
        x = round(atompos(2,tempi))-vol_size:round(atompos(2,tempi))+vol_size;
        y = round(atompos(1,tempi))-vol_size:round(atompos(1,tempi))+vol_size;
        z = round(atompos(3,tempi))-vol_size:round(atompos(3,tempi))+vol_size;
        atom_vol = atom_vol + volume(y,x,z).*filter;
        count = count + 1;
    end
end
atom_vol = atom_vol./count;

%% substact the traced atom
TraceVol = volume;
for tempi = 1:size(atompos,2)
    if ((round(atompos(2,tempi))-vol_size)>0) && ((round(atompos(1,tempi))-vol_size)>0) && ((round(atompos(3,tempi))-vol_size)>0) ...
       && (round(atompos(2,tempi))+vol_size<size(volume,2)) && (round(atompos(1,tempi))+vol_size<size(volume,1)) && (round(atompos(3,tempi))+vol_size<size(volume,3))
        x = round(atompos(2,tempi))-vol_size:round(atompos(2,tempi))+vol_size;
        y = round(atompos(1,tempi))-vol_size:round(atompos(1,tempi))+vol_size;
        z = round(atompos(3,tempi))-vol_size:round(atompos(3,tempi))+vol_size;
        TraceVol(y,x,z) = TraceVol(y,x,z) - atom_vol;
    end
end

%% convolution to generate new volume
% generate a 3D Gaussian filter to smooth the volume 
filter = atom_vol.*create3DGeneralizedGaussianFilter(vol_size,sigma_filter,2);filter = filter./sum(filter(:));
TraceVol = convn(TraceVol,filter,'same');
end


function filter = create3DGeneralizedGaussianFilter(size, sigma, beta)
    % size is the length of each dimension of the filter
    % sigma is the standard deviation of the Gaussian distribution

    % Create a grid of coordinates
    [x, y, z] = ndgrid(-size:size);

    % Calculate the Gaussian function
    filter = exp(-abs((x.^2 + y.^2 + z.^2) / (beta * sigma^2)).^beta);

    % Normalize the filter so that the sum of all elements is 1
    filter = filter ./ sum(filter(:));
end