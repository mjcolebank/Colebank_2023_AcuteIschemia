% Author: MJ Colebank
% Last edited: July 2023
%
% THIS CODE COMES WITH NO GUARANTEES
% Run a local sensitivity analysis of the triseg model parameterized for
% mouse hemodynamics
clear; clc; close all;
%%
for mouseID = 1:3
    NC     = 30;
    % [IC,pars] = get_pars_MI(mouseID,1.0); % Get the initial conditions and parameter values
    [IC,pars] = get_pars_NEW(mouseID,1.0); % Get the initial conditions and parameter values
    Tend   = NC.*pars.T;
    dt     = 1e-4;
    % pars.CV(9) = pars.CV(9)./2;

    if mouseID==1
        load GLS_5var_newSA_m1_5_24_23.mat q_hist par_ids parsnew
        qopt = q_hist{end};
        parsnew(par_ids) = exp(qopt);
        load opt_IC_m1.mat IC
        fname = 'sens_m1_final';
    elseif mouseID==2
         load GLS_5var_newSA_m2_5_24_23.mat q_hist par_ids parsnew
        qopt = q_hist{end};
        parsnew(par_ids) = exp(qopt);
        load opt_IC_m2.mat IC
        fname = 'sens_m2_final';
    else
         load GLS_5var_newSA_m3_5_25_23.mat q_hist par_ids parsnew 
        qopt = q_hist{end};
        parsnew(par_ids) = exp(qopt);
        load opt_IC_m3.mat IC
        fname = 'sens_m3_final';
    end
    pars.V     = parsnew(1:10);
    pars.SarcA = parsnew(11:25);
    pars.SarcV = parsnew(26:39);
    pars.CV    = parsnew(40:51);
    pars.Peri  = parsnew(52:53);
    pars.T     = parsnew(54);

    tspace = 0:dt:Tend;

    f = @(q) call_model_sens(q,IC,tspace,NC);
    
    STEP = 0.01;
    NORMALIZE = 0;
    LOG_PARS = 1;
%     par_ids = 1:53;
% par_ids=52
    [sens,nom_sol] = local_sensitivity_log(f,pars,STEP,NORMALIZE,par_ids,LOG_PARS);

    save(fname);
end