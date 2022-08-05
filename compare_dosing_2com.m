function compare_dosing_2com()
% COMPARE_DOSING-2COM() compares the plasma concentration time course 
% for a two compartment model. Treatment is administered as an infusion

%model parameters
p.CL   = 1.38E+01;   % central clearance
p.V1   = 1.48E+01;   % volume of distribution in central compartment 
p.Q    = 2.17E+00;   % inter-compartmental clearance
p.V2   = 4.23E+00;   % volume of distribution peripheral compartment
p.k    = p.CL/p.V1;  % rate constant of elimination              
p.k12  = p.Q/p.V1;   % rate constant from central to peripheral             
p.k21  = p.Q/p.V2;   % rate constant from peripheral to central           
p.regime=1;
p.endtime=168; %h

tspan = [0 p.endtime];     % max. time domain (h)
c0 = [0 0];         % Initial concentration in each compartment
tcourse=0:.01:tspan(2);


p.regime=2;
p.regime;
[t_vals_2,c_vals_2] = ode45(@derivatives, tspan, c0, [], p);
timedose2=r(tcourse,p);

p.regime=1;
[t_vals_1,c_vals_1] = ode45(@derivatives, tspan, c0, [], p);
timedose1=r(tcourse,p);


f = figure;
f.Position = [100 100 1050 400];
xlim_1=[0 tspan(2)];
xlim_1=[0 190];
subplot(3, 1, 1)
plot(tcourse,timedose2,'r')
hold on
plot(tcourse,timedose1,'b')
legend('regime2','regime1')
ylabel('\fontsize{13}Rate [mg/h]')
xlim(xlim_1)
subplot(3, 1, 2:3); 
plot(t_vals_1,c_vals_1(:,1),'-.b','DisplayName','regime 1:compartment 1 concentration')
hold on 
plot(t_vals_1,c_vals_1(:,2),'-b','DisplayName','regime 1:compartment 2 concentration')
plot(t_vals_2,c_vals_2(:,1),'-.r','DisplayName','regime 2:compartment 1 concentration')
plot(t_vals_2,c_vals_2(:,2),'-r','DisplayName','regime 2:compartment 2 concentration')
xlim(xlim_1)
ylim([0 6])
legend
xlabel('\fontsize{13}Time [hour]')
ylabel('\fontsize{13}Concentration [mg/L]')
%set(gca,"FontSize",10)

saveas(gcf,'concentration_profiles.png')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dcdt = derivatives(t, c, p)
dcdt = [r(t,p)/p.V1 - (p.k+p.k12)*c(1) + p.k21*p.V2/p.V1*c(2) 
        p.k12*p.V1/p.V2*c(1) - p.k21*c(2)                    ];          
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rt = r(t,p)
    switch p.regime
        case 1
        interval=72;
        duration=2;
        number_of_intervals=p.endtime/interval;
        sum_dosing=0;
        for counter0=0:1:number_of_intervals
            sum_dosing=sum_dosing+(counter0*interval<t & t<counter0*interval+duration);
            rt = 100/2*sum_dosing;
        end
        

        case 2
        interval=36;
        duration=2;
%         interval=24;
%         duration=4;
        number_of_intervals=p.endtime/interval;
%         while 0<t & t<2
%         number_of_intervals
%         end 

        sum_dosing=0;
        for counter1=0:1:number_of_intervals
            sum_dosing=sum_dosing+(counter1*interval<t & t<counter1*interval+duration);
            rt = 2*100/2*sum_dosing;
        end
        
        otherwise
        rt=0;    
    end 
end