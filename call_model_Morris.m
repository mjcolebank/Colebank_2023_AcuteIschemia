% Author: MJ Colebank
% Last edited: July 2023
%
% THIS CODE COMES WITH NO GUARANTEES
% Wrapper file to call the model for Morris screening
function [y_out,J] = call_model_Morris(pars_in,IC,tspace,par_ids,pars)
pars(par_ids) = pars_in;
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
    M(20,20) = 0; %DAE - ymd
    NC = 30;

    nt   = length(tspace);
    Tend = tspace(end);
    pars.dt = 1e-4;

    options=odeset('Mass',M,'RelTol',1e-12, 'AbsTol',1e-12,'MaxStep',1e-4);  %sets how accurate the ODE solver is
    tode = [0 0];
    ystorage = zeros(20,100*NC);
    tstorage = zeros(1,100*NC);
    for i=1:NC
        tode = [tode(2) i.*pars.T];
        ysol = ode15s(@DE_model_mass,tode,IC,options,pars);
        if isempty(ysol.y) || size(ysol.y,2)==1
            y_out = []; J=0;
            return
        elseif sum(isnan(ysol.y(:,end)))>0
            y_out = []; J=0;
            return
        end
        IC = ysol.y(:,end-1);
        tspace = linspace(ysol.x(1),ysol.x(end),100);
        ystorage(:,(i-1)*100+1:i*100) = deval(ysol,tspace);
        tstorage(:,(i-1)*100+1:i*100) = tspace;
    end
    if ysol.x(end)<Tend
        warning('ODE solver stopped prematurely. Recalculating.');
        %%%%%%%%%%%%%%%
        y_out = []; J=0;
        return
        %%%%%%%%%%%%%%
    end

else
    y_out  = [];
    J = 0;
    warning('Incompatible DAE conditions.');
    return;
end

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



