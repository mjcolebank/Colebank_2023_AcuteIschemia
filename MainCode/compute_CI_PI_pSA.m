% Author: MJ Colebank
% Last edited: July 2023
%
% THIS CODE COMES WITH NO GUARANTEES
% Compute the confidence intervals for the parameters based on the linear
% approximation to the error variance and parameter covariance matrix

function [upper_par, lower_par] = compute_CI_PI_pSA(par_opt,states,sens,data,W_err,V_NLS)
num_par = length(par_opt);
p_SA = states(1,:);
p_RV = states(2,:)';
V_RV = states(3,:)';
p_LV = states(4,:)';
V_LV = states(5,:)';

%% data
p_lv_data = data.p_lv_data;
p_rv_data = data.p_rv_data;

V_lv_data = data.V_lv_data;
V_rv_data = data.V_rv_data;

p_sa_data = data.p_sa_data;

n_cyc_LV = 1;
n_cyc_RV = 1;

n_data_LV = size(p_lv_data,1);
n_data_RV = size(p_rv_data,1);
n_total = 2.*n_data_RV.*n_cyc_RV + 3.*n_data_LV.*n_cyc_LV;

p_SA = interp1(linspace(0,1,length(p_SA)),p_SA,linspace(0,1,n_data_LV));
p_LV = interp1(linspace(0,1,length(p_LV)),p_LV,linspace(0,1,n_data_LV));
V_LV = interp1(linspace(0,1,length(V_LV)),V_LV,linspace(0,1,n_data_LV));
p_RV = interp1(linspace(0,1,length(p_RV)),p_RV,linspace(0,1,n_data_RV));
V_RV = interp1(linspace(0,1,length(V_RV)),V_RV,linspace(0,1,n_data_RV));


sens_p_SA = squeeze(sens(1,:,:));
sens_p_RV = squeeze(sens(2,:,:));
sens_p_LV = squeeze(sens(3,:,:));
sens_V_RV = squeeze(sens(4,:,:));
sens_V_LV = squeeze(sens(5,:,:));
sens_p_SA = interp1(linspace(0,1,length(sens_p_SA)),sens_p_SA,linspace(0,1,n_data_LV));
sens_p_RV = interp1(linspace(0,1,length(sens_p_RV)),sens_p_RV,linspace(0,1,n_data_RV));
sens_V_RV = interp1(linspace(0,1,length(sens_V_RV)),sens_V_RV,linspace(0,1,n_data_RV));
sens_p_LV = interp1(linspace(0,1,length(sens_p_LV)),sens_p_LV,linspace(0,1,n_data_LV));
sens_V_LV = interp1(linspace(0,1,length(sens_V_LV)),sens_V_LV,linspace(0,1,n_data_LV));

sens_full =[sens_p_RV; sens_V_RV; sens_p_LV; sens_V_LV; sens_p_SA];

% Dynamic residuals

r_LVP = (p_lv_data(:)-p_LV(:));
r_LVV = (V_lv_data(:)-V_LV(:));
r_RVP = (p_rv_data(:)-p_RV(:));
r_RVV = (V_rv_data(:)-V_RV(:));
r_PSA = (p_sa_data(:)-p_SA(:));

r_full = [r_RVP; r_RVV; r_LVP; r_LVV; r_PSA];%./sqrt(diag(W_err));

data_stack = [(data.p_rv_data); (data.V_rv_data); (data.p_lv_data); (data.V_lv_data); (data.p_sa_data)];
model_stack = [p_RV V_RV p_LV V_LV p_SA]';
%%

F_all    = sens_full'*inv(W_err)*sens_full;
% Get error variance estimator
s2 = diag(W_err);


if length(s2)==1
    s2 = s2.*ones(250,1);
end
invF_all_OG = inv(F_all);

% Covariance estimator
V_all    = s2(1).*invF_all_OG;%s2.*invF_all_OG;
if ~isempty(V_NLS)
    V_all = V_NLS;
end

SE2 = diag(V_all); % Variance for parameters
g = sens_full';%./exp(par_opt);       % sensitivity at each time point (NO WEIGHTS)

% Construct alpha/2 level confidence intervals
alpha = 0.975; % Two sided t-test

lower_par = par_opt - tinv(alpha,n_total-num_par).*sqrt(SE2);
upper_par = par_opt + tinv(alpha,n_total-num_par).*sqrt(SE2);


% Display the parameter confidence intervals and the +- values for paper
disp([par_opt tinv(alpha,n_total-num_par).*sqrt(SE2)])
%% Now construct model confidence and credible intervals
states_full = [p_RV; V_RV; p_LV; V_LV; p_SA];
CI_RVP = [p_RV; p_RV]';
CI_RVV = [V_RV; V_RV]';
CI_LVP = [p_LV; p_LV]';
CI_LVV = [V_LV; V_LV]';
CI_PSA = [p_SA; p_SA]';

PI_RVP = [p_RV; p_RV]';
PI_RVV = [V_RV; V_RV]';
PI_LVP = [p_LV; p_LV]';
PI_LVV = [V_LV; V_LV]';
PI_PSA = [p_SA; p_SA]';


n_RV = n_cyc_RV.*n_data_RV;
n_LV = n_cyc_LV.*n_data_LV;

t_val = tinv(alpha,n_total-num_par);
for i=1:n_total
    if i<=n_RV
        CI_RVP(i,1) = p_RV(i) - t_val.*sqrt(g(:,i)'*V_all*g(:,i));
        CI_RVP(i,2) = p_RV(i) + t_val.*sqrt(g(:,i)'*V_all*g(:,i));
        PI_RVP(i,1) = p_RV(i) - t_val.*sqrt(g(:,i)'*V_all*g(:,i)+s2(i));
        PI_RVP(i,2) = p_RV(i) + t_val.*sqrt(g(:,i)'*V_all*g(:,i)+s2(i));
    elseif i<=2.*n_RV
        CI_RVV(i-n_RV,1) = V_RV(i-n_RV) - t_val.*sqrt(g(:,i)'*V_all*g(:,i));
        CI_RVV(i-n_RV,2) = V_RV(i-n_RV) + t_val.*sqrt(g(:,i)'*V_all*g(:,i));
        PI_RVV(i-n_RV,1) = V_RV(i-n_RV) - t_val.*sqrt(g(:,i)'*V_all*g(:,i)+s2(i));
        PI_RVV(i-n_RV,2) = V_RV(i-n_RV) + t_val.*sqrt(g(:,i)'*V_all*g(:,i)+s2(i));
    elseif i<=2.*n_RV+n_LV
        CI_LVP(i-2.*n_RV,1) = p_LV(i-2.*n_RV) - t_val.*sqrt(g(:,i)'*V_all*g(:,i));
        CI_LVP(i-2.*n_RV,2) = p_LV(i-2.*n_RV) + t_val.*sqrt(g(:,i)'*V_all*g(:,i));
        PI_LVP(i-2.*n_RV,1) = p_LV(i-2.*n_RV) - t_val.*sqrt(g(:,i)'*V_all*g(:,i)+s2(i));
        PI_LVP(i-2.*n_RV,2) = p_LV(i-2.*n_RV) + t_val.*sqrt(g(:,i)'*V_all*g(:,i)+s2(i));
    elseif i<=2.*n_RV+2.*n_LV
       CI_LVV(i-2.*n_RV-n_LV,1) = V_LV(i-2.*n_RV-n_LV) - t_val.*sqrt(g(:,i)'*V_all*g(:,i));
        CI_LVV(i-2.*n_RV-n_LV,2) = V_LV(i-2.*n_RV-n_LV) + t_val.*sqrt(g(:,i)'*V_all*g(:,i));
        PI_LVV(i-2.*n_RV-n_LV,1) = V_LV(i-2.*n_RV-n_LV) - t_val.*sqrt(g(:,i)'*V_all*g(:,i)+s2(i));
        PI_LVV(i-2.*n_RV-n_LV,2) = V_LV(i-2.*n_RV-n_LV) + t_val.*sqrt(g(:,i)'*V_all*g(:,i)+s2(i));
    else
        CI_PSA(i-2.*n_RV-2.*n_LV,1) = p_SA(i-2.*n_RV-2.*n_LV) - t_val.*sqrt(g(:,i)'*V_all*g(:,i));
        CI_PSA(i-2.*n_RV-2.*n_LV,2) = p_SA(i-2.*n_RV-2.*n_LV) + t_val.*sqrt(g(:,i)'*V_all*g(:,i));
        PI_PSA(i-2.*n_RV-2.*n_LV,1) = p_SA(i-2.*n_RV-2.*n_LV) - t_val.*sqrt(g(:,i)'*V_all*g(:,i)+s2(i));
        PI_PSA(i-2.*n_RV-2.*n_LV,2) = p_SA(i-2.*n_RV-2.*n_LV) + t_val.*sqrt(g(:,i)'*V_all*g(:,i)+s2(i));
    end
end
%% Now plot prediction intervals
t_RV = linspace(0,data.T.*n_cyc_RV,n_RV)';
t_LV = linspace(0,data.T.*n_cyc_LV,n_LV)';
%%
figure(10);clf;
hold on;
fill([t_RV; flip(t_RV)],[PI_RVP(:,2); flip(PI_RVP(:,1))],[0.9 0.9 0.9],'EdgeAlpha',0);
fill([t_RV; flip(t_RV)],[CI_RVP(:,2); flip(CI_RVP(:,1))],[0.6 0.6 0.6],'EdgeAlpha',0);
plot(t_RV,p_RV,'r','LineWidth',3);
plot(t_RV,p_rv_data(:),'k','LineWidth',3)

figure(20);clf;
hold on;
fill([t_RV; flip(t_RV)],[PI_RVV(:,2); flip(PI_RVV(:,1))],[0.9 0.9 0.9],'EdgeAlpha',0);
fill([t_RV; flip(t_RV)],[CI_RVV(:,2); flip(CI_RVV(:,1))],[0.6 0.6 0.6],'EdgeAlpha',0);
plot(t_RV,V_RV,'r','LineWidth',3);
plot(t_RV,V_rv_data(:),'k','LineWidth',3)


figure(30);clf;
hold on;
fill([t_LV; flip(t_LV)],[PI_LVP(:,2); flip(PI_LVP(:,1))],[0.9 0.9 0.9],'EdgeAlpha',0);
fill([t_LV; flip(t_LV)],[CI_LVP(:,2); flip(CI_LVP(:,1))],[0.6 0.6 0.6],'EdgeAlpha',0);
plot(t_LV,p_LV,'r','LineWidth',3);
plot(t_LV,p_lv_data(:),'k','LineWidth',3)

figure(40);clf;
hold on;
fill([t_LV; flip(t_LV)],[PI_LVV(:,2); flip(PI_LVV(:,1))],[0.9 0.9 0.9],'EdgeAlpha',0);
fill([t_LV; flip(t_LV)],[CI_LVV(:,2); flip(CI_LVV(:,1))],[0.6 0.6 0.6],'EdgeAlpha',0);
plot(t_LV,V_LV,'r','LineWidth',3);
plot(t_LV,V_lv_data(:),'k','LineWidth',3)

figure(50);clf;
hold on;
fill([t_LV; flip(t_LV)],[PI_PSA(:,2); flip(PI_PSA(:,1))],[0.9 0.9 0.9],'EdgeAlpha',0);
fill([t_LV; flip(t_LV)],[CI_PSA(:,2); flip(CI_PSA(:,1))],[0.6 0.6 0.6],'EdgeAlpha',0);
plot(t_LV,p_SA,'r','LineWidth',3);
plot(t_LV,p_sa_data(:),'k','LineWidth',3)

end