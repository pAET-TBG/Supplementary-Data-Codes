function R1 = R1_objective_bayes(T, xdata, projections)
    persistent cacheKey cacheVal
    if isempty(cacheKey)
        cacheKey = strings(0);
        cacheVal = [];
    end

    p = [T.p11 T.p12 T.p13 T.p21 T.p22 T.p23];
    key = sprintf('%.6g_', p);  % rounding key

    hit = find(cacheKey == key, 1);
    if ~isempty(hit)
        R1 = cacheVal(hit);
        return
    end

    % ---- compute ----
    para = reshape(p,2,3);
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
        if ~isfinite(R1), R1 = 1e30; end
    catch
        R1 = 1e30;
    end

    cacheKey(end+1) = key;
    cacheVal(end+1,1) = R1;
end