function [cal_proj_opt,para_opt] = Cal_HBproj_bayes(para0, xdata, projections)

Z_arr	   = xdata.Z_arr;
Res        = xdata.Res;
halfWidth  = xdata.halfWidth;
% step_sz    = xdata.step_sz;

model      = xdata.model;
angles     = xdata.angles;
atom       = xdata.atoms;
num_atom   = numel(atom);
atom_type_num = numel(unique(atom));

[N1,N2,Num_pj] = size(projections);
%num_pj=size(angles,1);
% N_s = 2*halfWidth+1;


%para0 = [0.00047412, 0.00066645, 0.0011446;  33.774, 55.004, 62.278];
para0=para0;
lb = para0 * 0.8;%0.95;
ub = para0 * 1.2;%1.05;
% ===== define optimizable variables =====
%%
p = gcp('nocreate');
if isempty(p)
    parpool('local');  % or parpool(N)
end
%% ---- optimizable variables (explicit indexing, no reshape) ----
vars = [
    optimizableVariable('p11',[lb(1,1) ub(1,1)],'Type','real')
    optimizableVariable('p21',[lb(2,1) ub(2,1)],'Type','real')
];
%% ---- Bayesian optimization (parallel evaluations) ----
addpath src
results = bayesopt(@(T)objective_fun_par(T,xdata,projections), ...
                   vars, ...
                   'MaxObjectiveEvaluations',200, ...          % increase if needed
                   'NumSeedPoints',50, ...
                   'AcquisitionFunctionName','expected-improvement-plus', ...
                   'IsObjectiveDeterministic',true, ...
                   'UseParallel',true, ...                    % <<<<<< parallel
                   'ParallelMethod','clipped-model-prediction', ... % good default
                   'GPActiveSetSize',300, ...                 % optional speed-up
                   'PlotFcn',{@plotMinObjective,@plotObjectiveModel}, ...
                   'Verbose',1);

%% ---- best parameters ----
bestT = results.XAtMinObjective;
para_opt = [bestT.p11;
            bestT.p21];

R_opt = results.MinObjective;

disp('Best para =')
disp(para_opt)
disp('Best Rfactor (Objective) =')
disp(R_opt)

%% ---- final forward model ----
%%
[cal_proj_opt,~] = Cal_Bproj_2type(para_opt, xdata, projections);


for i=1:Num_pj
    pj = projections(:,:,i);
    resi_i=projections(:,:,i)-cal_proj_opt(:,:,i);
    Rarr(i) = sum(abs(resi_i(:)))/ sum(abs(pj(:)));

    pj_cal = cal_proj_opt(:,:,i);
    pj_cal = ImageNorm(pj_cal,sum(sum(pj)));
    resi_i_norm=pj-pj_cal;
    Rarr_norm(i) = sum(abs(resi_i_norm(:)))/ sum(abs(pj(:)));

    clear resi_i  pj  pj_cal resi_i
end

R=mean(Rarr)
R_norm=mean(Rarr_norm)
end