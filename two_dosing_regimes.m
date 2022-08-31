
function two_dosing_regimes()
% two_dosing_regimes compares the plasma concentration time course 
% for a two compartment model. Treatment is administered as an pill

%model parameters
% p.CL   = 0.01;     % central clearance
% p.V1   = 1;    % volume of distribution in central compartment 
% p.Q    = 0.15;   % inter-compartmental clearance
% p.V2   = 1;   % volume of distribution peripheral compartment
% p.ka= 0.7
% p.F=0.5

p.CL   = 0.3;     % central clearance
p.V1   = 1;    % volume of distribution in central compartment 
p.Q    = 0.15;   % inter-compartmental clearance
p.V2   = 1;   % volume of distribution peripheral compartment
p.ka= 0.7
p.F=0.5


p.k_elim    = p.CL/p.V1;  % rate constant of elimination              
p.k12  = p.Q/p.V1;   % rate constant from central to peripheral             
p.k21  = p.Q/p.V2;   % rate constant from peripheral to central   

p.min_effective_conc=500
p.max_tolerated_conc=700
p.regime=[];

p.maintence_dose=[];
p.load_dose=[];
p.interval=[];%%% prev value =24/3 

p.endtime=24*7*6; %h (6 weeks)

[p.interval, p.maintence_dose,max_dose] = dosing_chooser_blank( @derivatives,p.min_effective_conc, p.max_tolerated_conc,p)
p.load_dose=max_dose

for regime=[1 2]
    switch regime
        case 1
        p.regime=1;
        p.initial_dose=p.load_dose
        case 2
        p.regime=2;
        p.initial_dose=p.maintence_dose
        otherwise
            %
    end
    [t_vals_whole, c_vals_whole]=setandrunODE(p); 

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = figure;

f.Position = [100 100 1050 400];
unit_conv=24
xlim_1=[0 p.endtime/unit_conv]
%xlim_1=[0 72/unit_conv]
xtick_spacing=p.interval*3/unit_conv
plottime= linspace(0,p.endtime,100)

%xlim_1=[0 190];

%plot(t_vals_1/unit_conv,c_vals_1(:,1),'-.k','DisplayName','drug amount')
hold on 
plot(t_vals_1/unit_conv,c_vals_1(:,2),'-.b','DisplayName','regime 1:compartment 1 concentration')

plot(t_vals_1/unit_conv,c_vals_1(:,3),'-b','DisplayName','regime 1:compartment 2 concentration')
%plot(t_vals_2/unit_conv,c_vals_2(:,1),'-.k','DisplayName','drug amount')
plot(t_vals_2/unit_conv,c_vals_2(:,2),'-.r','DisplayName','regime 2:compartment 1 concentration')
plot(t_vals_2/unit_conv,c_vals_2(:,3),'-r','DisplayName','regime 2:compartment 2 concentration')
plot(plottime,p.min_effective_conc*ones(size(plottime)),'-c','DisplayName','min effective conc')
plot(plottime,p.max_tolerated_conc*ones(size(plottime)),'-c','DisplayName','max tolerated conc')

xlim(xlim_1);
xticks([xlim_1(1):xtick_spacing:xlim_1(2)])
ylim([0 p.max_tolerated_conc*1.2])
title( 'Regime 1 (Loading Dose) and Regime 2 (No Loading Dose)')
legend('Location','northeast')
xlabel('\fontsize{13}Time [day]')
ylabel('\fontsize{13}Concentration [mg/L]')
%set(gca,"FontSize",10)

saveas(gcf,'two_dosing_regimes_plot.png')

end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t_vals_whole, c_vals_whole]=setandrunODE(p)   
    %set up integration
    t_vals_whole=[];
    c_vals_whole=[];
    tspan = [0 p.interval];     % max. time domain (h)
    c0 =[p.initial_dose,0, 0]
    options = odeset('RelTol',1e-12, 'AbsTol',[1e-12 1e-12 1e-12]);
    %run stepwise integration
    while tspan(1)<p.endtime
        [t_vals,c_vals]= ode45(@derivatives, tspan, c0, options, p);
        tspan = [t_vals(end) t_vals(end)+p.interval]; 
        c0 = [c_vals(end,1)+p.maintence_dose, c_vals(end,2), c_vals(end,3)];
        t_vals_whole=vertcat(t_vals_whole,t_vals);
        c_vals_whole=vertcat(c_vals_whole,c_vals);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dcdt = derivatives(t, c, p)
a=c(1);
c1=c(2);
c2=c(3);

% dcdt = DERIVS_LinearAbs_Linear_Elim__rate_constants(a,c1,c2,params)

dcdt = [ - 1*p.ka*a, 
         p.F*p.ka*(a/p.V1) - p.k_elim*c1 - p.k12*c1*p.V1 + p.k21*c2*p.V2/p.V1,        
                                         +p.k12*c1*p.V1/p.V2 -  p.k21*c2*p.V2 ];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dose_interval, dose_amount,max_dose] = dosing_chooser_blank(derivativefunction, min_effective_conc, max_tolerated_conc,p)
    dose_interval=24/3
    dose_amount=2000
    max_dose=3000
end 
