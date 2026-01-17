%% The 3D and 2D Hanning window (Hanning filter) for 3D and 2D Fourier Transform
function KSpaceout = HanningFFT(array)
    [len,wid,hei]=size(array);
    HanningWindow=zeros(len,wid,hei);
    for i=1:len
        for j=1:wid
            for k=1:hei
                if hei > 1
                    HanningWindow(i,j,k)=0.125*(1-cos(2*pi*(i-1)/(len-1)))*(1-cos(2*pi*(j-1)/(wid-1)))*(1-cos(2*pi*(k-1)/(hei-1)));
                elseif hei == 1
                    HanningWindow(i,j,k)=0.125*(1-cos(2*pi*(i-1)/(len-1)))*(1-cos(2*pi*(j-1)/(wid-1)));
                end
            end
        end
    end
    array=array.*HanningWindow;
    KSpaceout = fftshift(fftn((ifftshift(array))));
end