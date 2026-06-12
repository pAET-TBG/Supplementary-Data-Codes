function R1 = objective_fun(T,xdata,projections)

    para = [T.p11 T.p12 T.p13;
            T.p21 T.p22 T.p23];

    try
        [cal_proj,~] = Cal_Bproj_2type_YXZ(para, xdata, projections);

        Num_pj = size(projections,3);
        Rarr = zeros(Num_pj,1);

        for i = 1:Num_pj
            pj = projections(:,:,i);
            resi_i = pj - cal_proj(:,:,i);
            Rarr(i) = sum(abs(resi_i(:))) / (sum(abs(pj(:))) + eps);
        end

        R1 = mean(Rarr);

        % 保护：bayesopt 不接受 NaN/Inf
        if ~isfinite(R1) || ~isscalar(R1)
            R1 = 1e30;
        end

    catch
        R1 = 1e30;
    end

%     % 可选：每一步打印（不影响 bayesopt 表格）
%     fprintf('R1=%.6g | [%.4g %.4g %.4g ; %.4g %.4g %.4g]\n', ...
%         R1, para(1,1), para(1,2), para(1,3), para(2,1), para(2,2), para(2,3));
end