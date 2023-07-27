% Author: MJ Colebank
% Last edited: July 2023
%
% THIS CODE COMES WITH NO GUARANTEES
% Wrapper file to call the model for optimization
% Returns the negative log-likelihood under the assumption of independent,
% mulitvariate Gaussian noise
function J = call_model_GLS(q,pars0,par_ids,data,IC,tspace,rnom,cov)
if ~isstruct(pars0)
    pars = pars0;
    pars(par_ids) = exp(q);
    pars = change_pars(pars);
else
    pars.V     = pars.V;%exp(pars.V);
    pars.SarcA = pars.SarcA;%exp(pars.SarcA);
    pars.SarcV = pars.SarcV;%exp(pars.SarcV);
    pars.CV    = pars.CV;%exp(pars.CV);
    pars.T     = pars.T;%exp(pars.T);
end
NC = 30;


%%
% NOTE: since we are changing parameters, the IC for the DAE need to be
% recalibrated, especially when changing sarcomere dynamics
%%
Vm_LV = - IC(8) - 0.5.*(pars.V(2)+pars.V(5))+IC(9);
Vm_RV =   IC(4) + 0.5.*(pars.V(4)+pars.V(5))+IC(9);
Vm = [Vm_LV Vm_RV IC(9)];
[xm,Am,Cm] = get_wallsegment(Vm,IC(20));
[ynew,Vm,Tx0,Ty0] = get_initial_conditions(0,xm,Vm,IC(20),IC,pars,[IC(8);IC(4)]);

if abs(Tx0)<1e-6 && abs(Ty0)<1e-6

    IC(9) = Vm(3);
    IC(20) = ynew;
    %% MASS MATRIX DAE APPROACH
    M = eye(20);
    M(9,9)   = 0; %DAE - Vsw
    M(20,20) = 0; %DAE - ym

    nt   = length(tspace);
    Tend = 0;
    pars.dt = 1e-4;

    options=odeset('InitialStep',1e-8,'Mass',M,'RelTol',1e-10, 'AbsTol',1e-10);  %sets how accurate the ODE solver is
     for i=1:NC
        Tstart = Tend;
        Tend = i.*pars.T;
        ysol = ode15s(@DE_model_mass,[Tstart Tend],IC,options,pars);
        IC  = ysol.y(:,end);
    end
    Bad_sol = 0;
    if ysol.x(end)<tspace(end)
        warning('ODE solver stopped prematurely. Recalculating.');
        Bad_sol=1;
    end
else
    Bad_sol=1;
end


% Load in the data + error variance
p_lv_data = data.p_lv_data;
p_rv_data = data.p_rv_data;
p_sa_data = data.p_sa_data;

V_lv_data = data.V_lv_data;
V_rv_data = data.V_rv_data;
nt_data = length(p_sa_data);



if Bad_sol==1
    
    J = rnom;
    return;
end

tstart = max(find(tspace>pars.T*(NC-1),1),find(tspace<=pars.T,1));
toutput = linspace(ysol.x(1),ysol.x(end),50);
ysols = deval(ysol,toutput);

out = get_model_results(toutput,ysols,pars);

p  = out.p./0.133322;%133.322;
V  = out.V.*1e3;

% Check for negative volumess
temp = V(1:8,:);
if any(temp(:)<0)
    
    J = rnom;
    return;
end


p_LV = p(8,:);
V_LV = V(8,:);
p_RV = p(4,:);
V_RV = V(4,:);
p_SA = p(1,:);

%%
% Dynamic residuals
r_LVP = (p_LV(:)-p_lv_data(:));
r_LVV = (V_LV(:)-V_lv_data(:));
r_RVP = (p_RV(:)-p_rv_data(:));
r_RVV = (V_RV(:)-V_rv_data(:));
r_SAP = (p_SA(:)-p_sa_data(:));


rfull = [r_RVP; r_RVV; r_LVP; r_LVV; r_SAP];

% Now return the negative log-likelihood with a penalty for the covariance
% scaling
log_det = log(det(cov));
if log_det<1e16
    J = log(det(cov))+rfull'*inv(cov)*rfull; % Neg-log likelihood
else
    J = rfull'*inv(cov)*rfull;
end
end

function pars = change_pars(pars)
temp.V     = pars(1:10);
temp.SarcA = pars(11:25);
temp.SarcV = pars(26:39);
temp.CV    = pars(40:51);
temp.Peri  = pars(52:53);
temp.T     = pars(54);
pars       = temp;
end