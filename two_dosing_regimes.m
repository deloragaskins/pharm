
function two_dosing_regimes()
% two_dosing_regimes compares the plasma concentration time course 
% for a two compartment model. Treatment is administered as an pill

%model parameters
p.CL   = 5;     % central clearance
p.V1   = 60;    % volume of distribution in central compartment 
p.Q    = 2.17;   % inter-compartmental clearance
p.V2   = 4.23;   % volume of distribution peripheral compartment
p.k    = p.CL/p.V1;  % rate constant of elimination              
p.k12  = p.Q/p.V1;   % rate constant from central to peripheral             
p.k21  = p.Q/p.V2;   % rate constant from peripheral to central   
p.ka= 15
p.regime=1;
p.endtime=24*7*6; %h
p.load_dose=1000
p.maintence_dose=100


for regime=[1 2]
    switch regime
        case 1
        p.regime=2;
        p.interval=24/3;
        p.initial_dose=p.load_dose
        case 2
        p.regime=2;
        p.interval=24/3;
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
xtick_spacing=p.interval*3/unit_conv
%xlim_1=[0 190];

plot(t_vals_1/unit_conv,c_vals_1(:,1),'-.b','DisplayName','regime 1:compartment 1 concentration')
hold on 
plot(t_vals_1/unit_conv,c_vals_1(:,2),'-b','DisplayName','regime 1:compartment 2 concentration')
plot(t_vals_2/unit_conv,c_vals_2(:,1),'-.r','DisplayName','regime 2:compartment 1 concentration')
plot(t_vals_2/unit_conv,c_vals_2(:,2),'-r','DisplayName','regime 2:compartment 2 concentration')
xlim(xlim_1);
xticks([xlim_1(1):xtick_spacing:xlim_1(2)])
ylim([0 p.load_dose*5/4])
title( 'Regime 1 (Loading Dose) and Regime 2 (No Loading Dose)')
legend('Location','southeast')
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
    c0 =[0, 0, p.initial_dose]
    options = odeset('RelTol',1e-12, 'AbsTol',[1e-12 1e-12 1e-12]);
    %run stepwise integration
    while tspan(1)<p.endtime
        [t_vals,c_vals]= ode45(@derivatives, tspan, c0, options, p);
        tspan = [t_vals(end) t_vals(end)+p.interval]; 
        c0 = [c_vals(end,1) c_vals(end,2), c_vals(end,3)+p.maintence_dose];
        t_vals_whole=vertcat(t_vals_whole,t_vals);
        c_vals_whole=vertcat(c_vals_whole,c_vals);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dcdt = derivatives(t, c, p)
dcdt = [ +p.ka*c(3)-p.k*c(1)+p.k12*c(1) + p.k21*p.V2/p.V1*c(2), 
                          -p.k21*c(2)+ p.k12*p.V1/p.V2*c(1),        
         -p.ka*c(3)                                                    ];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%