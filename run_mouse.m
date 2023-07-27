% Driver file for running the Triseg model used in Colebank et al. 2023
% Author: MJ Colebank
% Last edited: July 2023
%
% THIS CODE COMES WITH NO GUARANTEES
% Abbreviations: LA -> left atrium, LV -> left ventricle, RA -> right
% atrium, RV -> right ventricle, S -> septum,
% Am -> midwall area, Cm -> midwall curvature, ef -> wall strain
clear; clc; close all;
%% Parameters and ICs
% Setting opt=0 will run the model with nominal parameters; opt=1 will
% instead load the optimal parameters as documented in the manuscript
% Note that optimal parameters are log-scaled, and must be exponetiated.
opt=1;
for mouseID = 1:3
    % Get the initial conditions and parameter values
    [IC,pars,data] = get_pars_NEW(mouseID,1.0);
    if opt==1
        if mouseID==1
            load matfiles_to_load/opt_mouse1 q_hist par_ids parsnew pars0 q_sample
        elseif mouseID==2
            load matfiles_to_load/opt_mouse2 q_hist par_ids parsnew pars0 q_sample
        else
            load matfiles_to_load/opt_mouse3 q_hist par_ids parsnew pars0 q_sample
        end
        % Save the last iteration of parameter values and re-establish the
        % parameter structure for the model
        qopt = q_hist{end};
        parsnew(par_ids) = exp(qopt);

        pars.V     = parsnew(1:10);
        pars.SarcA = parsnew(11:25);
        pars.SarcV = parsnew(26:39);
        pars.CV    = parsnew(40:51);
        pars.Peri  = parsnew(52:53);
        pars.T     = parsnew(54);

        % Reguess DAE conditions; these should be based on the parameter
        % values provided to the model

        % To get initial values for the septal location states and the sarcomere
        % states, we need to solve for zero tension to begin with.
        Vm_LV = - IC(8) - 0.5.*(pars.V(2)+pars.V(5))+IC(9);
        Vm_RV =   IC(4) + 0.5.*(pars.V(4)+pars.V(5))+IC(9);
        Vm = [Vm_LV Vm_RV IC(9)];
        [xm,~,~] = get_wallsegment(Vm,IC(20));
        [ym0,Vm,dY,dV] = get_initial_conditions(0,xm,Vm,IC(20),IC,pars,IC([8 4]));
        [xm,Am,Cm] = get_wallsegment(Vm,ym0);

        IC(9)  = Vm(3);
        IC(20) = ym0;

    end
    %% Next, define the time vector
    num_cycles = 30;                 % Number of cycles to run the model for
    dt         = 1e-4;               % Time stepping for the output and plotting
    Tend       = num_cycles.*pars.T; % The end time point for the ODE/DAE solver
    tspace     = 0:dt:Tend;          % The time vector
    pars.dt    = dt;                 % Append the time stepping

    %% Now solve the model
    % To include DAE's, a mass matrix is defined that is the identity with
    % two diagonal elements being 0 reflecting the DAE components
    M = eye(20);  % Define the mass matrix
    M(9,9)   = 0; % DAE - Vsw
    M(20,20) = 0; % DAE - ym

    options=odeset('Mass',M,'RelTol',1e-10, 'AbsTol',1e-10);  %sets how accurate the ODE solver is

    %% Timing statements and model solve using ODE15s
    tic
    tode = [0 0];
    ystorage = zeros(20,100*num_cycles);
    tstorage = zeros(1,100*num_cycles);
    for i=1:num_cycles %Cycle through heart beats
        tode = [tode(2) i.*pars.T];
        ysol = ode15s(@DE_model_mass,tode,IC,options,pars);
        IC = ysol.y(:,end-1);
        tspace = linspace(ysol.x(1),ysol.x(end),100);
        ystorage(:,(i-1)*100+1:i*100) = deval(ysol,tspace);
        tstorage(:,(i-1)*100+1:i*100) = tspace;
    end
    toc

    %% Obtain all relevant outputs (volumes, pressures, sarcomere dynamics)
    n_cyc_plot = 2; %Number of cycles you want to plot
    tstart = find(tstorage>pars.T*(num_cycles-n_cyc_plot),1);
    toutput = tstorage(tstart:end);
    youtput = ystorage(:,tstart:end);

    % Get model outputs
    tic
    outputs = get_model_results(toutput,youtput,pars);
    toc


    %% Now plot outputs
    p      = outputs.p./0.133322;  % KPa  --> mmHg  (see get_model_results for output order)
    V      = outputs.V.*1e3;       % mL   --> muL   (see get_model_results for output order)
    q      = outputs.q.*1e3;       % mL/s --> muL/s (see get_model_results for output order)
    Gamma  = outputs.Sarc(1:5,:);  % Contractility (Gamma in the manuscript) (LA,LV,RA,RV,S)
    SL     = outputs.Sarc(6:10,:); % Sarcomere contractile length (Lsc)      (LA,LV,RA,RV,S)
    Am     = outputs.Am;           % Midwall area of cardiac chambers        (LA,LV,RA,RV,S)
    Cm     = outputs.Cm;           % Curvature of cardiac chambers           (LA,LV,RA,RV,S)
    ef     = outputs.ef;           % Strain of the cardiac chamber           (LA,LV,RA,RV,S)
    stress = outputs.stress;       % Wall stress of the cardiac chamber      (LA,LV,RA,RV,S)

    % Time vector for plotting purposes
    t      = toutput - toutput(1);
    toutput = toutput-toutput(1);

    %% Load data
    % Plot full time-series data (with "_stack" variables as well as the
    % mean signals that are used for model calibration
    load(strcat('matfiles_to_load/m',num2str(mouseID),'_timeseries_final'))
    load(strcat('matfiles_to_load/m',num2str(mouseID),'_mean_var'))
   

    p_LV = mu_PLV;
    V_LV = mu_VLV;
    p_RV = mu_PRV;
    V_RV = mu_VRV;
    p_SA = mu_PSA;

    %% Plotting
    figure(mouseID);
    subplot(2,2,1); hold on;
    plot(V(7,:),p(7,:),'LineWidth',3);
    title('LA'); set(gca,'FontSize',20); ylabel('Pressure (mmHg'); grid on;

    subplot(2,2,2); hold on;
    plot(V(3,:),p(3,:),'LineWidth',3);
    title('RA');set(gca,'FontSize',20); grid on;

    subplot(2,2,3); hold on;
    plot(VLV_stack(:),PLV_stack(:),'LineWidth',2,'Color',[0.5 0.5 0.5])
    plot(V_LV,p_LV,'k','LineWidth',2);
    plot(V(8,:),p(8,:),'LineWidth',3);
    ylabel('Pressure (mmHg'); xlabel('Volume (ml)');
    title('LV'); set(gca,'FontSize',20); grid on;

    subplot(2,2,4); hold on;
    plot(VRV_stack(:),PRV_stack(:),'LineWidth',2,'Color',[0.5 0.5 0.5])
    plot(V_RV,p_RV,'k','LineWidth',2);
    plot(V(4,:),p(4,:),'LineWidth',3);
    xlabel('Volume (ml)'); grid on;
    title('RV'); set(gca,'FontSize',20);
    %%
    figure(10+mouseID); hold on;
    plot(t(1:2:100),PSA_stack,'Color',[0.7 0.7 0.7]);
    plot(t(1:2:100),p(1,1:2:100),'c','LineWidth',3);
    ylabel('SA Pressure (mmHg)'); xlabel('Time (s)'); grid on;
    set(gca,'FontSize',20); grid on;

    %% Plot the time-series data
    figure(100+mouseID); clf;
    t = linspace(0,pars.T,length(p_LV));
    subplot(2,2,1); hold on;
    plot(linspace(0,pars.T.*2,2.*size(PLV_stack,1)),[PLV_stack; PLV_stack],'Color',[0.5 0.5 0.5],'LineWidth',1.5)
    plot([t t+t(end)],[p_LV; p_LV],'k','LineWidth',3);
    plot(toutput-toutput(1),p(8,:),'-.r','LineWidth',3);
    grid on;
    ylabel('LV pressure (mmHg)');
    set(gca,'FontSize',20);

    subplot(2,2,2); hold on;
    plot(linspace(0,pars.T.*2,2.*size(PLV_stack,1)),[VLV_stack; VLV_stack],'Color',[0.5 0.5 0.5],'LineWidth',1.5)
    plot([t t+t(end)],[V_LV; V_LV],'k','LineWidth',3);
    plot(toutput-toutput(1),V(8,:),'-.r','LineWidth',3);
    grid on;
    ylabel('LV Volume (\mu l)');
    set(gca,'FontSize',20);

    % RV
    subplot(2,2,3); hold on;
    plot(linspace(0,pars.T.*2,2.*size(PLV_stack,1)),[PRV_stack; PRV_stack],'Color',[0.5 0.5 0.5],'LineWidth',1.5)
    plot([t t+t(end)],[p_RV; p_RV],'k','LineWidth',3);
    plot(toutput-toutput(1),p(4,:),'-.r','LineWidth',3);
    grid on;
    ylabel('RV pressure (mmHg)');
    xlabel('Time (s)')
    set(gca,'FontSize',20);

    subplot(2,2,4); hold on;
    plot(linspace(0,pars.T.*2,2.*size(PLV_stack,1)),[VRV_stack; VRV_stack],'Color',[0.4 0.4 0.4],'LineWidth',1.5)
    plot([t t+t(end)],[V_RV; V_RV],'k','LineWidth',3);
    plot(toutput-toutput(1),V(4,:),'-.r','LineWidth',3);
    grid on;
    ylabel('RV Volume (\mu l)');
    xlabel('Time (s)')
    set(gca,'FontSize',20);

    %%
    figure(999);
    subplot(2,3,mouseID); hold on;
    plot(VLV_stack,PLV_stack,'Color',[0.6 0.6 0.6]);
    plot(V_LV,p_LV,'k','LineWidth',3);
    plot(V(8,:),p(8,:),'r','LineWidth',3);
    ylabel('Pressure (mmHg');
    xlabel('Volume (ml)');
    title('LV');
    set(gca,'FontSize',20);

    % RV
    subplot(2,3,3+mouseID); hold on;
    plot(VRV_stack,PRV_stack,'Color',[0.6 0.6 0.6]);
    plot(V_RV,p_RV,'k','LineWidth',3);
    plot(V(4,:),p(4,:),'b','LineWidth',3);
    ylabel('Pressure (mmHg');
    xlabel('Volume (ml)');
    title('RV');
    set(gca,'FontSize',20);

    %% Strain plotting
    long_strain = (SL-SL(:,1))./SL(:,1);
    green_strain = 0.5.*(SL.^2./SL(:,1).^2 - 1.0);
    circ_strain = (Am-Am(:,1))./Am(:,1);


    figure(1990);
    subplot(3,3,mouseID);
    plot(toutput,long_strain([2 4 5],:)','LineWidth',3)
    grid on; set(gca,'FontSize',20);
    ylim([-0.15 0.05])

    subplot(3,3,mouseID+3);
    plot(toutput,circ_strain([2 4 5],:)','LineWidth',3)
    grid on; set(gca,'FontSize',20);
    ylim([-0.20 0.05])

    subplot(3,3,mouseID+6);
    plot(toutput,green_strain([2 4 5],:)','LineWidth',3)
    grid on; set(gca,'FontSize',20);
    ylim([-0.15 0.05])

    %% Plot whole time series

    V_long = V(:,1:floor(end/2)+1);
    p_long = p(:,1:floor(end/2)+1);
    n_cyc_LV = size(PLV_stack,2);
    n_cyc_RV = size(PRV_stack,2);
    T_LV = pars.T.*n_cyc_LV;
    T_RV = pars.T.*n_cyc_RV;

    figure(30000+mouseID);clf;
    subplot(2,2,1); hold on;
    plot(linspace(0,T_LV,length(VLV_stack(:))),VLV_stack(:),'Color',[0.4 0.4 0.4],'LineWidth',2);
    plot(linspace(0,T_LV,size(V_long,2).*n_cyc_LV),repmat(V_long(8,:),1,n_cyc_LV),'r','LineWidth',2);

    subplot(2,2,2); hold on;
    plot(linspace(0,T_RV,size(V_long,2).*n_cyc_RV),repmat(V_long(4,:),1,n_cyc_RV),'r','LineWidth',2);
    plot(linspace(0,T_RV,length(VRV_stack(:))),VRV_stack(:),'Color',[0.4 0.4 0.4],'LineWidth',2);

    subplot(2,2,3); hold on;
    plot(linspace(0,T_LV,size(V_long,2).*n_cyc_LV),repmat(p_long(8,:),1,n_cyc_LV),'r','LineWidth',2);
    plot(linspace(0,T_LV,length(PLV_stack(:))),PLV_stack(:),'Color',[0.4 0.4 0.4],'LineWidth',2);

    subplot(2,2,4); hold on;
    plot(linspace(0,T_RV,size(V_long,2).*n_cyc_RV),repmat(p_long(4,:),1,n_cyc_RV),'r','LineWidth',2);
    plot(linspace(0,T_RV,length(PRV_stack(:))),PRV_stack(:),'Color',[0.4 0.4 0.4],'LineWidth',2);

end
