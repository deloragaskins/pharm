function pk_plus_hill_pd()
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
p.n0=1e9
p.kn_e=0
p.kn_a=4.5e-5
p.n_carry=1e12
p.kn_kill50=200
p.kn_killmax=0.01
p.n_factor=12




p.regime=2;
p.interval=24/3;
p.initial_dose=p.load_dose

[t_vals_whole, y_vals_whole]=setandrunODE(p); 
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = figure;

f.Position = [100 00 1000 400*2+40];
unit_conv=24
xlim_1=[0 p.endtime/unit_conv]
xtick_spacing=p.interval*3/unit_conv
%xlim_1=[0 190];

subplot(4, 1, 1:2)

plot(t_vals_whole/unit_conv,y_vals_whole(:,1),'-.b','DisplayName','compartment 1 concentration')
hold on 
% plot(t_vals_whole/unit_conv,y_vals_whole(:,2),'-b','DisplayName','compartment 2 concentration')
xlim(xlim_1);
xticks([xlim_1(1):xtick_spacing:xlim_1(2)])
ylim([0 p.load_dose*5/4])
legend('Location','southeast')
xlabel('\fontsize{13}Time [day]')
ylabel('\fontsize{13}Concentration [mg/L]')

title( 'Drug Concentration in compartment 1')

subplot(4, 1, 3:4)
plot(t_vals_whole/unit_conv,y_vals_whole(:,4),'-b','DisplayName','n')
ylim([0 p.n0*5/4])
xlim(xlim_1);
xticks([xlim_1(1):xtick_spacing:xlim_1(2)])


title( '___')
legend('Location','southeast')
xlabel('\fontsize{13}Time [day]')
ylabel('\fontsize{13} n')
%set(gca,"FontSize",10)

saveas(gcf,'pk_plus_hill_pd_plot.png')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t_vals_whole, y_vals_whole]=setandrunODE(p)   
    %set up integration
    t_vals_whole=[];
    y_vals_whole=[];
    tspan = [0 p.interval];     % max. time domain (h)
    y0 =[0, 0, p.initial_dose, p.n0]
    options = odeset('RelTol',1e-12, 'AbsTol',[1e-12 1e-12 1e-12 1e-12]);
    %run stepwise integration
    while tspan(1)<p.endtime
        [t_vals,y_vals]= ode45(@derivatives, tspan, y0, options, p);
        tspan = [t_vals(end) t_vals(end)+p.interval]; 
        y0 = [y_vals(end,1) y_vals(end,2), y_vals(end,3)+p.maintence_dose,y_vals(end,4)];
        t_vals_whole=vertcat(t_vals_whole,t_vals);
        y_vals_whole=vertcat(y_vals_whole,y_vals);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydt = derivatives(t, y, p)

kn_kill=HillEffect(p.n_factor*y(1),0,p.kn_killmax,p.kn_kill50,1);


dydt = [ +p.ka*y(3)-p.k*y(1)+p.k12*y(1) + p.k21*p.V2/p.V1*y(2), 
                          -p.k21*y(2)+ p.k12*p.V1/p.V2*y(1),        
         -p.ka*y(3),
         [(p.kn_a*log(p.n_carry/y(4)))-(p.kn_e+kn_kill)]*y(4)
                                                              ];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function E=HillEffect(c,E0,Emax,EC50,Hcoeff)
    E=E0+Emax.*c.^Hcoeff/(EC50.^Hcoeff+c.^Hcoeff);
end