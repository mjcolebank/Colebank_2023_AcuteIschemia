% Function that takes in the the sarcomere fiber strain and computes the
% total stress as a function of both active and pass stresses.

function [G_total,dGamma_dt,dLsc_dt] = sarcomere(t,eps_f,Lsc,Gamma,pars,AV)
active_scale = 1.0;
if AV==1
    Spars = pars.SarcA;
elseif AV==-1
    Spars = pars.SarcV;
    if length(pars.MI)>1
        active_scale = induce_MI(pars,t);%max(Spars(end-1)*0.1,Spars(end)); %take max of 30% active or the passive
    else
        active_scale = pars.MI;
    end
else
    Spars = pars.SarcV;
end
T     = pars.T;
tmod = mod(t,T);


Ls_ref      = Spars(1);
Ls_iso      = Spars(2);
V0          = Spars(3);
Lsc0        = Spars(4);
Gamma_rest  = Spars(5);

% Timing parameters are treated as parameters here: published code uses
% timing relative to some activation time
tauR  = Spars(6);
tauD  = Spars(7);
tauSC = Spars(8);

% Active and passive fiber stress values
sig_act    = Spars(9);
Ls_ref_pas = Spars(10);
sig_ECM    = Spars(11);
sig_Titin  = Spars(12);
k_ECM      = Spars(13);
k_Titin    = Spars(14);


% First, calculate the sarcomere length as a function of Sarcomere strain
Ls = Ls_ref.*exp(eps_f); %eq. B1

% Next, define the differential equation for the contractile element
% length, Lsc
dLsc_dt = ((Ls-Lsc)./Ls_iso - 1).*V0; %eq. B2

%% Active force
% Now, define the terms needed for the differential equation for sarcomere
% mechanical activiation, C

CL = tanh(4.0.*(Lsc-Lsc0).^2); %eq. B4

xF = min(8,tmod/tauR); % eq. B5b
F_rise = active_scale.*(0.02.*xF.^3 .*(8-xF).^2 .* exp(-xF)); %eq. B5

T_Lsc = tauSC.*(0.29+0.3.*Lsc); %eq. B6

% NOTE: This function has a discontinuity at t=T. Lumens actually uses an
% approximation (not included in the text;
exp_func = 1.0+exp((T_Lsc - tmod)./tauD); %in eq. B3

phi_rise  = CL.*F_rise./tauR;
phi_decay = ((Gamma_rest-Gamma)./exp_func)./tauD;
dGamma_dt = phi_rise + phi_decay; %eq. B3
%

% Calculate the active stress
G_act  = sig_act.*Gamma.*(Lsc-Lsc0).*(Ls-Lsc)./Ls_iso;

%% Passive force (titin and ECM)
% Now determine the passive stress, which depends on the contributions of
% titin and collagen to passive myocardial stiffness
% this component is based on the work from Vas Osta et al. 2020 and
% Walmsley et al. 2015

%Stretch in the passive components of the wall
Ls_pas = Ls./Ls_ref_pas; % Note: Ls = Ls_ref.*exp(eps_f)

% First, determine the stress due to ECM
G_ecm = max(sig_ECM.*(Ls_pas.^k_ECM - 1.0),0);

% Next, determine the stress due to titin, which also depends on active
% force
G_titin = max(sig_Titin.*(Ls_pas.^k_Titin - 1.0),0);

G_total = G_ecm + G_titin + G_act;

end