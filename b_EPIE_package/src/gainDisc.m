function disc = gainDisc(sizeP,radius, varargin)
%MASKDISC creates a disc shaped mask of given size and radius. 1 inside the
%disc and 0 outside.
%   radius: radius of the probe (aperture)
%   disc: a N1-by-N2 array of the disc centered at each specified
%
%Optional Inputs
%   sizeP: a 1-by-2 row vector, [N1, N2], OR a # specifying the rows and 
%   columns of the output. By default, sizeP == ceil(2.4*radius).
%
%Remarks
%   see test function
%
%   last updated 10/06/2016


% sizeP = ceil(2.4*radius);
if nargin >= 3 && ~isempty(varargin{1})
    sizeP = varargin{1};
end

if length(sizeP) == 1
    sizeP = [sizeP, sizeP];
end

pos = ceil((sizeP+1)/2);
N1 = sizeP(1);
N2 = sizeP(2);

[col, row] = meshgrid(1:N2,1:N1);
disc = (((row - pos(1)).^2 + (col - pos(2)).^2) <= radius^2);


end
