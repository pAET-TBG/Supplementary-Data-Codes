function R1 = objective_fun_par(T,xdata,projections)

    para = [T.p11 T.p12;
            T.p21 T.p22];

    try
        [cal_proj,~] = Cal_Bproj_2type(para, xdata, projections);

        Num_pj = size(projections,3);
        Rarr = zeros(Num_pj,1);

        for i = 1:Num_pj
            pj = projections(:,:,i);
            resi_i = pj - cal_proj(:,:,i);
            Rarr(i) = sum(abs(resi_i(:))) / (sum(abs(pj(:))) + eps);
        end

        R1 = mean(Rarr);

        % bayesopt cannot accept NaN/Inf
        if ~isfinite(R1) || ~isscalar(R1)
            R1 = 1e30;
        end

    catch
        R1 = 1e30;
    end

    % (optional) print each evaluation:
    % fprintf('R1=%.6g | [%.4g %.4g %.4g ; %.4g %.4g %.4g]\n', ...
    %    R1, para(1,1), para(1,2), para(1,3), para(2,1), para(2,2), para(2,3));
end