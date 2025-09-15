function peak_point = fit_Guassian_function(ydata,xdata,flag)
% normalize the xdata and ydata
xdata_org = xdata;
xdata     = xdata  - min(xdata(:));
xdata     = xdata  ./ max(xdata(:));
ydata     = ydata' - min(ydata(:));
ydata     = ydata  ./ max(ydata(:));

% set the parameter
opts = optimset('Display','off');
fun  = @(p,xdata) abs(p(1))*exp(-((xdata-p(2))/p(3)).^2);
p0   = [1,0.5,length(xdata)/2^5];  % the initial guess for the fitting
lb   = [0,   0, 0];
ub   = [inf, 1, inf];

% fit with Gaussian
[p,~] = lsqcurvefit(fun,p0,xdata,double(ydata),lb,ub,opts);

% de-normalization
peak_point = p(2);
peak_point = peak_point .* max(xdata_org(:)-min(xdata_org(:)));
peak_point = peak_point + min(xdata_org(:));

% visualize the fitting result
if flag == 1
    yfit = fun(p,xdata);
    figure();plot(xdata,ydata);hold on;
    plot(xdata,yfit);hold off
end

end