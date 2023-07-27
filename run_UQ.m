% Author: MJ Colebank
% Last edited: July 2023
%
% THIS CODE COMES WITH NO GUARANTEES
% Construct confidence and prediction intervals
clear; clc; close all;
%% Start with parameters
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
%% Load files
which_mouse = 1;
folder_name = 'matfiles_to_load/';
load(strcat(folder_name,'opt_mouse',num2str(which_mouse)),'q_hist','Hessian_hist','V_hist')
load(strcat(folder_name,'opt_mouse',num2str(which_mouse),'_outputs'),'outputs','data','pars','p','V')
res_sens = load(strcat(folder_name,'sens_res_m',num2str(which_mouse),'.mat'),'sens');
load(strcat(folder_name,'sens_output_m',num2str(which_mouse)),'sens','par_ids')


qopt = q_hist{end};
W_err = V_hist{end};
S_res = squeeze(res_sens.sens);
V_cov = inv(S_res'*inv(W_err)*S_res);

%%
if ~isempty(V_cov)
    figure(100);
    c = zeros(11,11);
    name_cell = {};
    for i=1:11
        for j=1:11
            c(i,j) = V_cov(i,j)./sqrt(V_cov(i,i)*V_cov(j,j));
        end
        name_cell{end+1} = Names(par_ids(i));
    end
    H = figure(99); 
    h=heatmap(triu(c));
    h.XDisplayLabels= name_cell;
    h.YDisplayLabels= name_cell;
    h.CellLabelFormat = '%0.2f';
    colormap(jet);
    set(gca,'FontSize',20);
    H.Renderer = "painters";
    H.Position = [4.1860e+02   1.2306e+03   8.6160e+02   5.7280e+02];
%     print(strcat(fname,'_corr'),'-depsc')
end
%% Look at mean signal 
data2.n_cyc_LV = 1;
data2.n_cyc_RV = 1;
data2.p_lv_data = mean(data.p_lv_data,2);
data2.V_lv_data = mean(data.V_lv_data,2);
data2.p_rv_data = mean(data.p_rv_data,2);
data2.V_rv_data = mean(data.V_rv_data,2);
data2.p_sa_data = mean(data.p_sa_data,2);


states = [p(1,1:99); p(4,1:99); V(4,1:99); p(8,1:99); V(8,1:99)];
sens = sens([1 4 8 13 17],:,:);
data2.T = pars.T;
[CI_upp,CI_low] = compute_CI_PI_pSA(qopt,states,sens,data2,W_err,V_cov);
par_CI = exp([CI_low qopt CI_upp]);
[exp(CI_low)./exp(qopt) exp(CI_upp)./exp(qopt)]


%%
% return;
h = figure(10); axis tight; axis tight; grid on; set(gca,'FontSize',20);
h.Position = [1.0000e+00   4.1000e+01   1.2800e+03   6.0733e+02];

h = figure(20); axis tight;axis tight; grid on; set(gca,'FontSize',20);
h.Position = [1.0000e+00   4.1000e+01   1.2800e+03   6.0733e+02];

h = figure(30); axis tight; axis tight; grid on; set(gca,'FontSize',20);
h.Position = [1.0000e+00   4.1000e+01   1.2800e+03   6.0733e+02];

h = figure(40); axis tight; grid on; set(gca,'FontSize',20);
h.Position = [1.0000e+00   4.1000e+01   1.2800e+03   6.0733e+02];

h = figure(50); axis tight; grid on; set(gca,'FontSize',20);
h.Position = [1.0000e+00   4.1000e+01   1.2800e+03   6.0733e+02];

