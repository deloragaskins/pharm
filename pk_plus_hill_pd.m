function pk_plus_hill_pd()
% 
addpath(genpath('./Toolbox')) %required for derivatives function
derivativefunction=@DERIVS_2comp_LinearAbs_Linear_Elim__rate_constants_hill;

%model parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%
p.n0=1e9;
deriv_params.n_carry=1e12;
deriv_params.kn_e=0;
deriv_params.kn_a=4.5e-5;
deriv_params.kn_kill50=200;
deriv_params.kn_killmax=0.01;
deriv_params.n_factor=12;
%%%%%%%%%%%%
%p.n_factor
%p.kn_killmax
%p.kn_kill50
%kn_a
%p.kn_e

p.load_dose=1000;
p.maintence_dose=100;
p.interval=24/3;

p.endtime=24*7*6; %h
%run time
%course%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[t_vals_whole, y_vals_whole]=run_dosing_course__two_compartment_PD(derivativefunction, p.endtime,p.interval,p.load_dose,p.maintence_dose,p.n0,deriv_params);     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = figure;
f.Position = [100 00 1000 400*2+40];
unit_conv=24;
xlim_1=[0 p.endtime/unit_conv];
xtick_spacing=p.interval*3/unit_conv;
%xlim_1=[0 190];

subplot(4, 1, 1:2)

plot(t_vals_whole/unit_conv,y_vals_whole(:,1),'-.b','DisplayName','compartment 1 concentration')
hold on 
% plot(t_vals_whole/unit_conv,y_vals_whole(:,2),'-b','DisplayName','compartment 2 concentration')
xlim(xlim_1);
xticks(xlim_1(1):xtick_spacing:xlim_1(2))
ylim([0 p.load_dose*5/4])
legend('Location','southeast')
xlabel('\fontsize{13}Time [day]')
ylabel('\fontsize{13}Concentration [mg/L]')

title( 'Drug Concentration in compartment 1')

subplot(4, 1, 3:4)
plot(t_vals_whole/unit_conv,y_vals_whole(:,4),'-b','DisplayName','n')
ylim([0 p.n0*5/4])
xlim(xlim_1);
xticks(xlim_1(1):xtick_spacing:xlim_1(2))


title( '___')
legend('Location','southeast')
xlabel('\fontsize{13}Time [day]')
ylabel('\fontsize{13} n')
%set(gca,"FontSize",10)

saveas(gcf,'pk_plus_hill_pd_plot.png')
end

