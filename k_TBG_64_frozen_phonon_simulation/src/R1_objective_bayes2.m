function T_out = R1_objective_bayes2(T, xdata, projections)

    para = [T.p11 T.p12 T.p13;
            T.p21 T.p22 T.p23];

    [cal_proj,~] = Cal_Bproj_2type_YXZ(para, xdata, projections);

    Num_pj = size(projections,3);
    Rarr = zeros(Num_pj,1);
    for i = 1:Num_pj
        pj = projections(:,:,i);
        resi_i = pj - cal_proj(:,:,i);
        Rarr(i) = sum(abs(resi_i(:))) / (sum(abs(pj(:))) + eps);
    end

    R1 = mean(Rarr);
    if ~isfinite(R1), R1 = 1e30; end

    % 让表格列名显示 Rfactor（可选）
    T_out = table(R1,'VariableNames',{'Rfactor'});
end