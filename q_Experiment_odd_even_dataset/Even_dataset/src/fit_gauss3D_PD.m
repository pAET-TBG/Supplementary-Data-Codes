
function [x, resnorm,residual] = fit_gauss3D_PD(init_guess, xdata, ydata, varargin)

%fit_gauss3D_PD(init_guess, xdata, ydata, fixed, lowerbound,upperboud)
% Y. Yang, UCLA, 2015
% This function is just after slight modification from
% original version by M. Bartels, UCLA, 2014
% To make the Gaussian always positively defined

%Fit 3D gaussian and allow fixe parameters as well as bounds

%set standard values for optional parameters
fixed=zeros(size(init_guess));
lowerbound= zeros(size(init_guess)) - Inf;
upperbound= zeros(size(init_guess)) + Inf;

%check optional parameters
if nargin>3
   fixed=cell2mat(varargin(1));
end

if nargin>4
   lowerbound=cell2mat(varargin(2));
end

if nargin>5
   upperbound=cell2mat(varargin(3));
end

%convert to double
init_guess    = double(init_guess);
ydata = double(ydata);

%store all initial parameters
init_guess_all = init_guess;

%only fit the variable parameters
init_guess = init_guess_all(fixed == 0);
lowerbound = lowerbound(fixed == 0);
upperbound = upperbound(fixed == 0);


opt=optimset('TolFun',1e-12);
opt=optimset(opt,'Display','off');

%do the fit
[x, resnorm,residual] = lsqcurvefit(@Fhelp, init_guess, xdata, ydata, lowerbound, upperbound, opt);

%add fixed parameters to the final result again
init_guess_all(fixed == 0) = x;
x = init_guess_all;

    
    function y = Fhelp(x, xdata)
        %helper function to deal with fixed parameters
        
        %merge with fixed parameters
        x_current = x;
        x = init_guess_all;
        x(fixed == 0) = x_current;
        
        %calculate the acutal function
        y=calc_gauss3D_PD(x,xdata);

    end

end

function y = calc_gauss3D_PD(x, xdata)
%Calulate 3D gaussian f(v) = exp(-v'*A*v)
%based on code by R.Xu, UCLA, 2014

    [L,M,N] = size(xdata.x); Num = L*M*N;

    v = [reshape(xdata.y - x(3),1,Num); reshape(xdata.x - x(4),1,Num); reshape(xdata.z - x(5),1,Num)];

    vector1 = [1 0 0];
    rotmat1 = MatrixQuaternionRot(vector1,x(9));
    
    vector2 = [0 1 0];
    rotmat2 = MatrixQuaternionRot(vector2,x(10));

    vector3 = [0 0 1];
    rotmat3 = MatrixQuaternionRot(vector3,x(11));

    rotMAT =  rotmat3*rotmat2*rotmat1;

    D = [1/x(6)  0   0;
         0  1/x(7)  0;
         0    0  1/x(8)];

    A = rotMAT' * D * rotMAT; 
    y = x(2)*reshape(exp(-dot(v,A*v,1)),L,M,N) + x(1);    
    
end
function dd = MatrixQuaternionRot(vector,theta)

% theta = theta*pi/180;
vector = vector/sqrt(dot(vector,vector));
w = cos(theta/2); x = -sin(theta/2)*vector(1); y = -sin(theta/2)*vector(2); z = -sin(theta/2)*vector(3);
RotM = [1-2*y^2-2*z^2 2*x*y+2*w*z 2*x*z-2*w*y;
      2*x*y-2*w*z 1-2*x^2-2*z^2 2*y*z+2*w*x;
      2*x*z+2*w*y 2*y*z-2*w*x 1-2*x^2-2*y^2;];

dd = RotM;
end