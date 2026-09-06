function [x, resnorm, residual] = fit_gauss1D_PD(init_guess, xdata, ydata, varargin)
% fit_gauss1D_PD(init_guess, xdata, ydata, fixed, lowerbound, upperbound)
%
% Fit a positively defined 1D Gaussian:
%
%   y = x(2) * exp(-(xdata-x(3)).^2 / x(4)) + x(1)
%
% Parameters:
%   x(1) = background
%   x(2) = amplitude
%   x(3) = center
%   x(4) = Gaussian width parameter
%
% Example:
%   init_guess = [0, 10, 25, 4];
%
% Optional:
%   fixed      = [0 0 0 0];      % 1 = fixed, 0 = fitted
%   lowerbound = [-Inf 0 -Inf 0];
%   upperbound = [ Inf Inf Inf Inf];


%% Default values
fixed      = zeros(size(init_guess));
lowerbound = zeros(size(init_guess)) - Inf;
upperbound = zeros(size(init_guess)) + Inf;


%% Optional inputs
if nargin > 3
    fixed = cell2mat(varargin(1));
end

if nargin > 4
    lowerbound = cell2mat(varargin(2));
end

if nargin > 5
    upperbound = cell2mat(varargin(3));
end


%% Convert to double
init_guess = double(init_guess);
xdata      = double(xdata);
ydata      = double(ydata);


%% Store all parameters
init_guess_all = init_guess;


%% Only fit variable parameters
init_guess = init_guess_all(fixed == 0);
lowerbound = lowerbound(fixed == 0);
upperbound = upperbound(fixed == 0);


%% Optimization settings
opt = optimset('TolFun', 1e-12);
opt = optimset(opt, 'Display', 'off');


%% Fit
[x, resnorm, residual] = lsqcurvefit( ...
    @Fhelp, ...
    init_guess, ...
    xdata, ...
    ydata, ...
    lowerbound, ...
    upperbound, ...
    opt);


%% Put fixed parameters back
init_guess_all(fixed == 0) = x;
x = init_guess_all;


%% Helper function
    function y = Fhelp(x_fit, xdata)

        % Start with all original parameters
        x_current = init_guess_all;

        % Replace fitted parameters
        x_current(fixed == 0) = x_fit;

        % Calculate Gaussian
        y = calc_gauss1D_PD(x_current, xdata);

    end

end


function y = calc_gauss1D_PD(x, xdata)
% Calculate 1D Gaussian
%
% x(1) = background
% x(2) = amplitude
% x(3) = center
% x(4) = width parameter

y = x(2) .* exp(-(xdata - x(3)).^2 ./ x(4)) + x(1);

end