%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|       Author: CXS group, PSI                                          |
%*-----------------------------------------------------------------------*
% You may use this code with the following provisions:
%
% If the code is fully or partially redistributed, or rewritten in another
%   computing language this notice should be included in the redistribution.
%
% If this code, or subfunctions or parts of it, is used for research in a 
%   publication or if it is fully or partially rewritten for another 
%   computing language the authors and institution should be acknowledged 
%   in written form in the publication: “Data processing was carried out 
%   using the “cSAXS matlab package” developed by the CXS group,
%   Paul Scherrer Institut, Switzerland.” 
%   Variations on the latter text can be incorporated upon discussion with 
%   the CXS group if needed to more specifically reflect the use of the package 
%   for the published work.
%
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%   
% This code and subroutines are part of a continuous development, they 
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its 
%    proper use and the correctness of the results.

function img = imshift_fft(img, x,y, apply_fft, split)
    % IMSHIFT_FFT  will apply subpixel shift that can be different for each frame. If apply_fft == false, then images will be assumed to be in fourier space 

    if nargin  < 3
        y = x(:,2);
        x = x(:,1);
    end
    if nargin < 4
        apply_fft = true;
    end
    if nargin < 5 
        split = 1; %% split for low GPU memory 
    end
    
    if all(x==0) && all(y==0)
        return
    end
    real_img = isreal(img);
    Np = size(img);
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     [nr,nc]=size(img);
%     Nr = ifftshift([-fix(nr/2):ceil(nr/2)-1]);
%     Nc = ifftshift([-fix(nc/2):ceil(nc/2)-1]);
%     [Nc,Nr] = meshgrid(Nc,Nr);
%     img2 = ifft2(fft2(img).*exp(2i*pi*(y*Nr/nr+x*Nc/nc)));
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if apply_fft
        if split == 1
            img = fft2(img);
        else
            img = fft2_partial(img);
        end
    end

    
    
    xgrid = (fftshift((0:Np(2)-1)/Np(2))-0.5);
    X = reshape((x(:)*xgrid)',1,Np(2),[]);
    X =  exp((-2i*pi)*X);
    img = bsxfun(@times, img,X);
    ygrid = (fftshift((0:Np(1)-1)/Np(1))-0.5);
    Y = reshape((y(:)*ygrid)',Np(1),1,[]);
    Y =  exp((-2i*pi)*Y);
    img = bsxfun(@times, img,Y);
        
    if apply_fft
        if split == 1
            img = ifft2(img);
        else
            img = ifft2_partial(img);
        end
    end
    if real_img
        img = real(img);
    end
    
  
end
