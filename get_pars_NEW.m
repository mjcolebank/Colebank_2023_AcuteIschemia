% Function that provides the nominal parameter values for the
% Triseg/circulation model.
% Author: MJ Colebank
% Last edited: 3/10/2022
%
% THIS CODE COMES WITH NO GUARANTEES
% Abbreviations: LA -> left atrium, LV -> left ventricle, RA -> right
% atrium, RV -> right ventricle, S -> septum, Am -> midwall area, Cm ->
% midwall curvature

function [IC,pars,data] = get_pars_NEW(mouseID,T_fac) %make variable input in future versions
% Set the nominal parameter values
if mouseID==1
    load matfiles_to_load/m1_mean_var.mat

    data.p_lv_data = mu_PLV;
    data.V_lv_data = mu_VLV;
    data.p_rv_data = mu_PRV;
    data.V_rv_data = mu_VRV;
    data.p_sa_data = mu_PSA;

    % Include the variances
    data.p_lv_var = var_PLV;
    data.p_rv_var = var_PRV;
    data.V_lv_var = var_VLV;
    data.V_rv_var = var_VRV;
    data.p_sa_var = var_PSA;


    T = 0.1034;

elseif mouseID==2
    load matfiles_to_load/m2_mean_var.mat
    data.p_lv_data = mu_PLV;
    data.V_lv_data = mu_VLV;
    data.p_rv_data = mu_PRV;
    data.V_rv_data = mu_VRV;
    data.p_sa_data = mu_PSA;

    % Include the variances
    data.p_lv_var = var_PLV;
    data.p_rv_var = var_PRV;
    data.V_lv_var = var_VLV;
    data.V_rv_var = var_VRV;
    data.p_sa_var = var_PSA;

    T = 0.1211;
else

    load matfiles_to_load/m3_mean_var.mat
    data.p_lv_data = mu_PLV;
    data.V_lv_data = mu_VLV;
    data.p_rv_data = mu_PRV;
    data.V_rv_data = mu_VRV;
    data.p_sa_data = mu_PSA;

    % Include the variances
    data.p_lv_var = var_PLV;
    data.p_rv_var = var_PRV;
    data.V_lv_var = var_VLV;
    data.V_rv_var = var_VRV;
    data.p_sa_var = var_PSA;

    T = 0.1129;

end
T = T*T_fac;
p_LV = mu_PLV;
V_LV = mu_VLV;
p_RV = mu_PRV;
V_RV = mu_VRV;
p_SA = mu_PSA;


%% Get the average heart rate

SV = max(V_LV)-min(V_LV);
SV = max(SV,(max(V_RV)-min(V_RV)));
SV = SV./1000;
HR = 60./T;
CO_microl_min = SV*HR;
CO_ml = CO_microl_min/60;

if mouseID==1
    Vtot = 84.7*29;
elseif mouseID==2
    Vtot = 84.7*32.4;
else
    Vtot = 84.7.*30.1;
end
Vtot = Vtot./1000;
%% Define conversion factors
p_conv = 0.133322;% mmHg -> Kpa

%% First, CV model parameters
% Most are from Boron and Boulapep, Tables 19-3 (vascular) and 22-3 (heart)

p_lv_sys  = max(p_LV)*p_conv;
p_lv_dias = min(p_LV)*p_conv;

p_rv_sys  = max(p_RV)*p_conv;
p_rv_dias = min(p_RV)*p_conv;

p_sa_sys  = 0.99*p_lv_sys;
p_sa_dias = min(p_SA).*p_conv;%0.66*p_lv_sys;%min(p_SA).*p_conv;%
p_sa_mean = (p_sa_sys+2.*p_sa_dias)./3;

p_pa_sys  = 0.99*p_rv_sys;
p_pa_dias = 0.32*p_rv_sys;%0.27*p_rv_sys;
p_pa_mean = (p_pa_sys+2.*p_pa_dias)./3;

p_sys_cap = 0.26*p_sa_mean;
p_pul_cap = 0.70*p_pa_mean;

p_sv_mean =  0.26*p_sa_mean;
p_pv_mean =  0.15*p_pa_mean;

p_la_dias = 0.45*p_pv_mean;

p_ra_dias = 0.45*p_sv_mean;


%% Now, initialize the volume estimates
% Using Lumens 2009 and Boron pg. 878
Vsa_tot   = 0.06*Vtot;%0.06*Vtot;   % mL
Vsv_tot   = 0.76*Vtot;   % mL
Vra_tot   = 0.010*Vtot;  % mL
Vrv_tot   = 0.026*Vtot;  % mL
Vpa_tot   = 0.027*Vtot;  % mL
Vpv_tot   = 0.072*Vtot;  % mL
Vla_tot   = 0.010*Vtot;  % mL
Vlv_tot   = 0.026*Vtot;  % mL
Vsw_tot   = 0.009*Vtot; % mL

Vm_SW = Vsw_tot;  % mL
ym    = 0.262957; % Initial guess for midwall junction (updated later)

%% Calculate stressed volumes based on Beneken and DeWitt 1966

Vsa = Vsa_tot.*0.27;
Vsv = Vsv_tot.*0.075;
Vra = Vra_tot.*1.0;
Vrv = Vrv_tot.*1.0;
Vpa = Vpa_tot.*0.58;
Vpv = Vpv_tot.*0.11;
Vla = Vla_tot.*1.0;
Vlv = Vlv_tot.*1.0;
Vsw = Vsw_tot.*1.0;

%% Next, define TriSeg parameters

% From Tim's data again!
if mouseID==1
    Am_ref_LA = 0.15;
    Am_ref_LV = 0.50;
    Am_ref_RA = 0.15;
    Am_ref_RV = 0.55;
    Am_ref_SW = 0.25;

    V_LAwall = 0.00285;
    V_LVwall = 0.05628;
    V_RAwall = 0.00199;
    V_RVwall = 0.01994;
    V_SWwall = 0.02841;
elseif mouseID==2
    Am_ref_LA = 0.15;
    Am_ref_LV = 0.70;
    Am_ref_RA = 0.15;
    Am_ref_RV = 0.75;
    Am_ref_SW = 0.25;

    V_LAwall = 0.00218;
    V_LVwall = 0.07673;
    V_RAwall = 0.00389;
    V_RVwall = 0.03713;
    V_SWwall = 0.03837;
else
    Am_ref_LA = 0.15;
    Am_ref_LV = 0.55;
    Am_ref_RA = 0.15;
    Am_ref_RV = 0.45;
    Am_ref_SW = 0.25;

    V_LAwall = 0.00389;
    V_LVwall = 0.05508;
    V_RAwall = 0.00299;
    V_RVwall = 0.01671;
    V_SWwall = 0.02754;

end
% Pericardial Volume
V0heart = (Vla+Vlv+Vra+Vrv).*0.9;


%% Parameters for the sarcomere model

% Atria
Ls_ref     = 2.0;         % micrometer
Ls_iso     = 0.04;        % micrometer
vmax       = 24;          % micrometer per second
Lsc0       = 1.51;        % micrometer
Gam_rest   = 0.02;        % dimensionless
tauR       = 0.0375*T;    % seconds
tauD       = 0.005*T;     % seconds
tauSC      = 0.15*T;      % seconds
sig_act    = 35;          % KPa
Ls_ref_pas = 1.8;         % micrometer
sig_ECM    = 0.08;        % Kpa
sig_Titin  = 0.25;        % Kpa
k_ECM      = 10;          % dimensionless
k_Titin    = 6;           % dimensionless

t_atr_offset = 0.18.*T;   % seconds

Spars_A = [Ls_ref; Ls_iso; vmax; Lsc0; Gam_rest; ...
    tauR; tauD; tauSC; sig_act; ...
    Ls_ref_pas; sig_ECM; sig_Titin; k_ECM; k_Titin; ...
    t_atr_offset];


% Ventricles
Ls_ref     = 2.0;         % micrometer
Ls_iso     = 0.04;        % micrometer
vmax       = 12;          % micrometer per second
Lsc0       = 1.51;        % micrometer
Gam_rest   = 0.02;        % dimensionless
tauR       = 0.009;       % seconds
tauD       = 0.009;       % seconds
tauSC      = 0.038;       % seconds
sig_act    = 75;          % KPa
Ls_ref_pas = 1.8;         % micrometer
sig_ECM    = 0.08;        % Kpa
sig_Titin  = 0.25;        % Kpa
k_ECM      = 10;          % dimensionless
k_Titin    = 6;           % dimensionless

Spars_V = [Ls_ref; Ls_iso; vmax; Lsc0; Gam_rest; ...
    tauR; tauD; tauSC; sig_act; ...
    Ls_ref_pas; sig_ECM; sig_Titin; k_ECM; k_Titin];

%% Pericardium
s = 10;
Peri = [V0heart; s];

%% Nominal resistance values
Ra_val = 0.1.*p_conv/CO_ml;              % KPa s / muL, Aortic Valve Resistance
Rm_val = 0.75.*p_conv/CO_ml;             % KPa s / muL, Mitral Valve Resistance
Rp_val = 0.1.*p_conv/CO_ml;              % KPa s / muL, Pulmonic Valve Resistance
Rt_val = 0.75.*p_conv/CO_ml;             % KPa s / muL, Tricuspid Valve Resistance
Rvc    = (p_sv_mean - p_ra_dias)./CO_ml; % KPa s / muL, Vena Cava Resistance
Rpv    = (p_pv_mean - p_la_dias)./CO_ml; % KPa s / muL, Pulmonary venous Resistance
Rs     = (p_sa_sys-p_sys_cap)./CO_ml;    % KPa s / muL, Systemic vascular Resistance
Rp     = (p_pa_sys-p_pul_cap)./CO_ml;    % KPa s / muL, Pulmonary vascular Resistance
Csa    = Vsa./(p_sa_sys);                % muL / KPa,   Systemic artery Compliance
Csv    = Vsv./p_sv_mean;                 % muL / KPa,   Systemic venous Compliance
Cpa    = Vpa./p_pa_sys;                  % muL / KPa,   Pulmonary artery Compliance
Cpv    = Vpv./p_pv_mean;                 % muL / KPa,   Pulmonary venous Compliance


% Set parameter values
Vpars = [V_LAwall; V_LVwall; V_RAwall; V_RVwall; V_SWwall; ...
    Am_ref_LA; Am_ref_LV; Am_ref_RA; Am_ref_RV; Am_ref_SW];


CVpars = [Ra_val; Rm_val; Rp_val; Rt_val; Rvc; Rpv; Rs; Rp; ...
    Csa; Csv; Cpa; Cpv];

if any(CVpars<0)
    error('Negative Parameters');
end

pars.V     = Vpars;
pars.SarcA = Spars_A;
pars.SarcV = Spars_V;
pars.CV    = CVpars;
pars.Peri  = Peri;
pars.T     = T;

%%
% To get initial values for the septal location states and the sarcomere
% states, we need to solve for zero tension to begin with.
Vm_LV = - Vlv - 0.5.*(V_LVwall+V_SWwall)+Vsw;
Vm_RV =   Vrv + 0.5.*(V_RVwall+V_SWwall)+Vsw;
Vm = [Vm_LV Vm_RV Vm_SW];
[xm,~,~] = get_wallsegment(Vm,ym);

Gam_la0 = Gam_rest;
Gam_lv0 = Gam_rest;
Gam_ra0 = Gam_rest;
Gam_rv0 = Gam_rest;
Gam_sw0 = Gam_rest;

Lla0 = Ls_ref;
Llv0 = Ls_ref;
Lra0 = Ls_ref;
Lrv0 = Ls_ref;
Lsw0 = Ls_ref;

ym0 = ym;

% Initial guess of initial conditions
IC = [Vsa; Vsv; Vra; Vrv; Vpa; Vpv; Vla; Vlv; Vsw;...
    Gam_la0; Gam_lv0; Gam_ra0; Gam_rv0; Gam_sw0; ...
    Lla0; Llv0; Lra0; Lrv0; Lsw0; ...
    ym0];

% Update the septal midwall volume and junction point by solving a root
% finding problem
% Vm = Vm([2 4 5]);
[ym0,Vm,~,~] = get_initial_conditions(0,xm,Vm,ym,IC,pars,[Vlv;Vrv]);
IC(9)  = Vm(3);
IC(20) = ym0;

end
