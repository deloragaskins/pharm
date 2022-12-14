
function two_dosing_regimes()
% two_dosing_regimes compares the plasma concentration time course 
% for a two compartment model. Treatment is administered as an pill


addpath(genpath('./Toolbox')) %required for dosing chooser and derivatives function
derivativefunction=@DERIVS_2comp_LinearAbs_Linear_Elim__rate_constants;

%model parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % p.CL   = 5;     % central clearance
% % % % % p.V1   = 60;    % volume of distribution in central compartment 
% % % % % p.Q    = 2.17;   % inter-compartmental clearance
% % % % % p.V2   = 4.23;   % volume of distribution peripheral compartment
% % % % % p.k    = p.CL/p.V1;  % rate constant of elimination              
% % % % % p.k12  = p.Q/p.V1;   % rate constant from central to peripheral             
% % % % % p.k21  = p.Q/p.V2;   % rate constant from peripheral to central   
% % % % % p.ka= 15
% % % % % p.regime=1;
% % % % % p.endtime=24*7*6; %h
% % % % % p.load_dose=1000
% % % % % p.maintence_dose=100

deriv_params.CL   = .693/24;     % central clearance
deriv_params.V1   = 1;    % volume of distribution in central compartment 
deriv_params.Q    = 0.15;   % inter-compartmental clearance
deriv_params.V2   = 1;   % volume of distribution peripheral compartment
deriv_params.ka= 0.7;
deriv_params.F=0.5;
%%%%%%%
deriv_params.k_elim    = deriv_params.CL/deriv_params.V1;  % rate constant of elimination              
deriv_params.k12  = deriv_params.Q/deriv_params.V1;   % rate constant from central to peripheral             
deriv_params.k21  = deriv_params.Q/deriv_params.V2;   % rate constant from peripheral to central   
%%%%%%%
p.min_effective_conc=500;
p.max_tolerated_conc=700;
p.regime=[];

p.maintence_dose=[];
p.load_dose=[];
p.interval=[];%%% prev value =24/3 

p.endtime=24*7*6; %h (6 weeks)

%obtain concentration vs time for 2 regimes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[p.interval, p.maintence_dose,max_dose] = dosing_chooser( derivativefunction,p.min_effective_conc, p.max_tolerated_conc,deriv_params);
p.load_dose=max_dose;

for regime=[1 2]
    switch regime
        case 1
        p.regime=1;
        p.initial_dose=p.load_dose;
        case 2
        p.regime=2;
        p.initial_dose=p.maintence_dose;
        otherwise
            %
    end
    [t_vals_whole, c_vals_whole]=run_dosing_course(derivativefunction,p.endtime,p.interval,p.initial_dose,p.maintence_dose, deriv_params);

    switch regime
        case 1
        t_vals_1= t_vals_whole;
        c_vals_1= c_vals_whole;
        case 2
        p.regime=2;
        t_vals_2= t_vals_whole;
        c_vals_2= c_vals_whole;
        otherwise      %
    end
    
end

%plot!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = figure;

f.Position = [100 100 1050 400];
unit_conv=24;
xlim_1=[0 p.endtime/unit_conv];
%xlim_1=[0 72/unit_conv]
%xlim_1=[0 190];
xtick_spacing=p.interval*3/unit_conv;
plottime= linspace(0,p.endtime,100);

%%%%
%plot(t_vals_1/unit_conv,c_vals_1(:,1),'-.k','DisplayName','drug amount')
hold on 
plot(t_vals_1/unit_conv,c_vals_1(:,2),'-.b','DisplayName','regime 1:compartment 1 concentration')
plot(t_vals_1/unit_conv,c_vals_1(:,3),'-b','DisplayName','regime 1:compartment 2 concentration')
%%%%
%plot(t_vals_2/unit_conv,c_vals_2(:,1),'-.k','DisplayName','drug amount')
plot(t_vals_2/unit_conv,c_vals_2(:,2),'-.r','DisplayName','regime 2:compartment 1 concentration')
plot(t_vals_2/unit_conv,c_vals_2(:,3),'-r','DisplayName','regime 2:compartment 2 concentration')
%%%%
plot(plottime,p.min_effective_conc*ones(size(plottime)),'-c','DisplayName','min effective conc')
plot(plottime,p.max_tolerated_conc*ones(size(plottime)),'-c','DisplayName','max tolerated conc')

xlim(xlim_1);
xticks(xlim_1(1):xtick_spacing:xlim_1(2))
ylim([0 p.max_tolerated_conc*1.2])

title( 'Regime 1 (Loading Dose) and Regime 2 (No Loading Dose)')
legend('Location','southeast')
xlabel('\fontsize{13}Time [day]')
ylabel('\fontsize{13}Concentration [mg/L]')
%set(gca,"FontSize",10)

saveas(gcf,'two_dosing_regimes_plot.png')

end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




