% Author: MJ Colebank
% Last edited: July 2023
%
% THIS CODE COMES WITH NO GUARANTEES
% Run an optimization on the TriSeg model with different noise variance
% Using the log-likelihood with MVN noise
clear; clc; close all;
Names = {'VW_{la}','VW_{lv}','VW_{ra}','VW_{rv}','VW_{s}',...
    'Amref_{la}','Amref_{lv}','Amref_{ra}','Amref_{rv}','Amref_{s}',...
    'Lsref_{a}','Lsiso_{a}','vmax_{a}','Lsc0_{a}','Crest_{a}',...
    'tauR_{a}','tauD_{a}','tauSC_{a}','sigact_{a}', ...
    'Lsref-pas_{a}','sig_{a,ECM}','sig_{a,Titin}','k_{a,ECM}','k_{a_Titin}',...
    'toffset',...
    'Lsref_{v}','Lsiso_{v}','vmax_{v}','Lsc0_{v}','Crest_{v}',...
    'tauR_{v}','tauD_{v}','tauSC_{v}','sigact_{v}', ...
    'Lsref-pas_{v}','sig_{v,ECM}','sig_{v,Titin}','k_{v,ECM}','k_{v_Titin}',...
    'Ra_val','Rm_val','Rp_val','Rt_val','Rvc','Rpv',...
    'Rs','Rp','Csa','Csv','Cpa','Cpv',...
    'V0','s'};
%%
mouseID = 1;
[IC,pars,data] = get_pars_NEW(mouseID,1.0); % Get the initial conditions and parameter values
% Select which parameters to infer
par_ids = [7 9 10 31 32 33 34 46 47 49 52];

n_par = length(par_ids);
pars0 = struct2cell(pars);
pars0 = cell2mat(pars0);
q0 = pars0(par_ids);

% do everything with a log-scale
upper_bounds = q0.*5;
lower_bounds = q0.*0.05;

q0 = log(q0);
upper_bounds = log(upper_bounds);
lower_bounds = log(lower_bounds);


%% Next, define the time vector
NC     = 30;
dt     = 5e-4;
Tend   = NC.*pars.T;
tspace = 0:dt:Tend;
pars.dt = dt;

%%




% Initial estimate of covariance for OLS to start
V0 = eye(250,250); % Could use data

Jnom = 1e8;%call_model_GLS(q0,pars0,par_ids,data,IC,tspace,0,V0);
F  = @(q) call_model_GLS(q,pars0,par_ids,data,IC,tspace,Jnom,V0);
res0 = get_res_GLS(q0,pars0,par_ids,data,IC,tspace,Jnom,V0);
q_tol = 1e-8;
curr_diff = 1;
V = V0;
% Now loop until convergence
q_hist = {};
V_hist = {};
Hessian_hist = {};
res_hist = {};
data_stack = [data.p_rv_data data.V_rv_data data.p_lv_data data.V_lv_data data.p_sa_data];

q_hist{end+1} = q0;
res_hist{end+1} = res0;
V_hist{end+1} = V0;
Hessian_hist{end+1} = 0;


% If using parallel
use_parallel = 0;
options=optimoptions('fminunc','Display','iter','Algorithm','quasi-newton',...
    'MaxFunctionEvaluations',80000,'MaxIterations',80000,'UseParallel',use_parallel==1);
if use_parallel==1
    parpool(6);
end

%%
% Try to use multiple initial conditions
n_q0_init = 1;
q0_orig = q0;
q_sample = q0_orig;
counter = 1;
qold = q0;
while curr_diff>q_tol
    F  = @(q) call_model_GLS(q,pars0,par_ids,data,IC,tspace,Jnom,V);

    [qopt,FVALupdate,EXITFLAG,OUTPUT,GRAD,HESSIAN] = fminunc(F,qold,options);

    curr_diff = sum((exp(qold)-exp(qopt)).^2);
    qold = qopt;
    parsnew = pars0;
    parsnew(par_ids) = exp(qopt);

    % Save things
    q_hist{end+1} = qopt;
    V_hist{end+1} = V;
    Hessian_hist{end+1} = HESSIAN;
    res_hist{end+1} = get_res_GLS(qopt,pars0,par_ids,data,IC,tspace,Jnom,V);
    res_test = get_res_GLS(qopt,pars0,par_ids,data,IC,tspace,Jnom,V0);

    % Obtain the diagonal variance estimate for each measuremenet source

    V(1:50,1:50)       = eye(50,50)*sum((res_test(1:50).^2))./50;
    V(51:100,51:100)   = eye(50,50)*sum((res_test(51:100).^2))./50;
    V(101:150,101:150) = eye(50,50)*sum((res_test(101:150).^2))./50;
    V(151:200,151:200) = eye(50,50)*sum((res_test(151:200).^2))./50;
    V(201:250,201:250) = eye(50,50)*sum((res_test(201:250).^2))./50;

end

save('GLS_results');
