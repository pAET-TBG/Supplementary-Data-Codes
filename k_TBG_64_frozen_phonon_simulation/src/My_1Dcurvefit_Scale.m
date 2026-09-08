% Hanfeng Zhong, UCLA
% 2024.11.15
function [p, fminres, fitresult, eflag] = My_1Dcurvefit_Scale(refdata, testdata,init_guess)
    opts = optimset('Display','off');
    % fitting for x axis projection
    fun = @(p,xdata) xdata*p(1);
    [p,fminres,~,eflag] = lsqcurvefit(fun,init_guess,testdata,refdata,[],[],opts);
    fitresult = testdata*p(1);
end