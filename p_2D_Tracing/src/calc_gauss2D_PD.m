function y = calc_gauss2D_PD(x,xdata)

% x(1): the offset for the 2D Gaussian
% x(2): the maximum number for the 2D Gaussian (the actual maximum is equal to x(1) + x(2))
% x(3): the center position for y-axis
% x(4): the center position for x-axis
% x(5): the sigma for x-x
% x(6): the sigma for y-y
% x(7): the rotation angle in 2D plane
% xdata: the coordinate for the 2D Gaussian

% Calculate 2D gaussian f(v) = exp(-v'*D*v)
% based on code by R.Xu, UCLA, 2014
[L, M] = size(xdata.x);
Num = L*M;
v = [reshape(xdata.y - x(3),1,Num);
     reshape(xdata.x - x(4),1,Num)];

D = [1/x(5) 0;
     0  1/x(6)];

theta = x(7);
rotMAT = [cos(theta), -sin(theta); sin(theta), cos(theta)];
A = rotMAT' * D * rotMAT; 

y = x(2)*reshape(exp(-dot(v,A*v,1)),L,M) + x(1);  

end