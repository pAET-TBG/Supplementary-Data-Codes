function para0 = calculate_B_factor(para0, xdata, model, projections)
opt = optimset('TolFun', 1e-12, 'TolX', 1e-8, 'MaxIter', 1000, 'Display', 'iter');
% subtract a constant background
projections(projections>0) = projections(projections>0) - min(projections(projections>0));
x0 = para0;                                                                % the parameter of the function perpared for optimization
x0(1,:)=x0(1,:)/x0(1,1);                                                   % normalize the para0(1,1) para0(1,2) para0(1,3)
xdata.model = model;                                                       % the positions of each atoms perpared for optimization
xdata.model_ori = model;                                                   % the original position of each atoms will not be changed during this iteration
xdata.projections=projections;
[para0, ~,~] = lsqcurvefit(@Cal_Bproj_2type2, x0, xdata, projections, [0 0; 0 0], [1 1 ;20 20], opt);
fprintf('H1 = %.03f, H2 = %.03f\n B1 = %.03f, B2 = %.03f\n', para0(1),para0(3),para0(2),para0(4));

end