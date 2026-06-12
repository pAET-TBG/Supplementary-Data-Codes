function [para_opt, R_opt, cal_proj_opt, out] = fit_para0_minR(para0, xdata, projections)
% Fit para0 (2x3) to minimize R1 = mean_i sum|pj-cal| / sum|pj|
%
% Inputs:
%   para0        - 2x3 initial parameters
%   xdata        - your xdata
%   projections  - Ny x Nx x Num_pj
%
% Outputs:
%   para_opt        - 2x3 optimized parameters
%   R_opt           - minimized R1
%   cal_proj_opt    - Cal_Bproj_2type_YXZ(para_opt, ...) output
%   out             - struct with optimizer info

    % ---- initial vector (6x1) ----
    p0 = para0(:);

    % ---- objective ----
    obj = @(p) R1_objective(p, xdata, projections);

    % ---- options ----
    opts = optimset('Display','iter', ...
                    'MaxIter', 50, ...
                    'MaxFunEvals', 3000, ...
                    'TolX', 1e-8, ...
                    'TolFun', 1e-10);

    % ---- run ----
    [p_opt, R_opt, exitflag, output] = fminsearch(obj, p0, opts);

    % ---- reshape back ----
    para_opt = reshape(p_opt, size(para0));
%para_opt=p_opt;
    % ---- final prediction ----
    [cal_proj_opt, ~] = Cal_Bproj_2type_YXZ(para_opt, xdata, projections);

    out.exitflag = exitflag;
    out.output   = output;
end


function R1 = R1_objective(p_vec, xdata, projections)
    % reshape to 2x3
    para = reshape(p_vec, 2, 3);

        [cal_proj, ~] = Cal_Bproj_2type_YXZ(para, xdata, projections);

        Num_pj = size(projections, 3);
        Rarr = zeros(Num_pj, 1);

        for i=1:Num_pj
            pj = projections(:,:,i);
            resi_i=projections(:,:,i)-cal_proj(:,:,i);
            Rarr(i) = sum(abs(resi_i(:)))/ sum(abs(pj(:)));
            clear resi_i  pj  pj_cal resi_i
        end

        R1 = mean(Rarr);
end