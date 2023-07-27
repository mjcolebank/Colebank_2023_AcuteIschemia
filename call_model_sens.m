% Author: MJ Colebank
% Last edited: July 2023
%
% THIS CODE COMES WITH NO GUARANTEES
% Wrapper file to call the model for local sensitivity analysis
function [y_out,J] = call_model_sens(pars,IC,tspace,NC)
if ~isstruct(pars)
    pars = change_pars(pars);
else
    pars.V     = pars.V;
    pars.SarcA = pars.SarcA;
    pars.SarcV = pars.SarcV;
    pars.CV    = pars.CV;
    pars.T     = pars.T;
end
%%
% NOTE: since we are changing parameters, the IC for the DAE need to be
% recalibrated, especially when changing sarcomere dynamics
%%
Vm_LV = - IC(8) - 0.5.*(pars.V(2)+pars.V(5))+IC(9);
Vm_RV =   IC(4) + 0.5.*(pars.V(4)+pars.V(5))+IC(9);
Vm = [Vm_LV Vm_RV IC(9)];
[xm,Am,Cm] = get_wallsegment(Vm,IC(20));
[ynew,Vm,~,~] = get_initial_conditions(0,xm,Vm,IC(20),IC,pars,[IC(8);IC(4)]);
IC(9) = Vm(3);
IC(20) = ynew;
%% MASS MATRIX DAE APPROACH
M = eye(20);
M(9,9)   = 0; %DAE - Vsw
M(20,20) = 0; %DAE - ym

nt   = length(tspace);
Tend = tspace(end);
pars.dt = 1e-4;

options=odeset('Mass',M,'RelTol',1e-8, 'AbsTol',1e-8);  %sets how accurate the ODE solver is
 tic
    tode = [0 0];
    ystorage = zeros(20,100*NC);
    tstorage = zeros(1,100*NC);
    for i=1:NC
        tode = [tode(2) i.*pars.T];
        ysol = ode15s(@DE_model_mass,tode,IC,options,pars);
        IC = ysol.y(:,end-1);
        tspace = linspace(ysol.x(1),ysol.x(end),100);
        ystorage(:,(i-1)*100+1:i*100) = deval(ysol,tspace);
        tstorage(:,(i-1)*100+1:i*100) = tspace;
    end
    toc


ysols = deval(ysol,tspace);

tstart = find(tstorage>pars.T*(NC-1),1);
toutput = tstorage(tstart:end);
youtput = ystorage(:,tstart:end);
out = get_model_results(toutput,youtput,pars);
p = out.p./0.133322;%133.322;
V = out.V.*1e3;
q = out.q.*1e3;
Gamsarc = out.Sarc(1:5,:);
Lsarc = out.Sarc(6:10,:);
Am = out.Am;
Cm = out.Cm;
ef = out.ef;
stress = out.stress;

y_out = [p; V; q; Gamsarc; Lsarc; Am; Cm; ef; stress];
J = 0;

end


function pars = change_pars(pars)
temp.V     = pars(1:10)';
temp.SarcA = pars(11:25)';
temp.SarcV = pars(26:39)';
temp.CV    = pars(40:51)';
temp.Peri  = pars(52:53)';
temp.T     = pars(54);
pars       = temp;
end

