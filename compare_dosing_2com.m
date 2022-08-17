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


for regime=[1 2]
    switch regime
        case 1
        p.regime=2;
        p.interval=48;
        p.duration=2;
        p.Rate_0 = 50*2;
        case 2
        p.regime=2;
        p.interval=24;
        p.duration=2;
        p.Rate_0 = 50;
        otherwise
            %
    end
    [t_vals_whole, c_vals_whole ,r_vals]=setandrunODE(p); 

    switch regime
        case 1
        t_vals_1= t_vals_whole;
        c_vals_1= c_vals_whole;
        r_vals_1=r_vals;
        case 2
        p.regime=2;
        t_vals_2= t_vals_whole;
        c_vals_2= c_vals_whole;
        r_vals_2=r_vals;
        otherwise      %
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = figure;
f.Position = [100 100 1050 400];
%xlim_1=[0 p.endtime];
xlim_1=[0 190];
subplot(3, 1, 1)
plot(t_vals_1,r_vals_1,'b','LineWidth',1,'DisplayName','regime1' )
hold on
plot(t_vals_2,r_vals_2,'r','DisplayName','regime2')
ylabel('\fontsize{13}Rate [mg/h]')
xlim(xlim_1);
subplot(3, 1, 2:3); 
plot(t_vals_1,c_vals_1(:,1),'-.b','DisplayName','regime 1:compartment 1 concentration')
hold on 
plot(t_vals_1,c_vals_1(:,2),'-b','DisplayName','regime 1:compartment 2 concentration')
plot(t_vals_2,c_vals_2(:,1),'-.r','DisplayName','regime 2:compartment 1 concentration')
plot(t_vals_2,c_vals_2(:,2),'-r','DisplayName','regime 2:compartment 2 concentration')
xlim(xlim_1);
ylim([0 6])
legend
xlabel('\fontsize{13}Time [hour]')
ylabel('\fontsize{13}Concentration [mg/L]')
%set(gca,"FontSize",10)

saveas(gcf,'compare_dosing_2com_plot.png')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t_vals_whole, c_vals_whole, r_vals]=setandrunODE(p)   
    %set up integration
    t_vals_whole=[];
    c_vals_whole=[];
    r_vals=[];
    ie=[1 2];
    
    p.Rate=p.Rate_0;
    stepsize=0.01;
    tspan = [0 : stepsize:p.endtime];     % max. time domain (h)
    c0 = [0 0];         % Initial concentration in each compartment
    options = odeset('Events',@dose_events,'RelTol',1e-12, 'AbsTol',[1e-12 1e-12]);
    %run stepwise integration
    while ~isempty(ie)
        if tspan(1)== p.endtime
            "oopsies" %bug catching
            %this happens because the while loop is not catching the last
            %loops empty ie
            break
        end
        [t_vals,c_vals,te,ye,ie]= ode15s(@derivatives, tspan, c0, options, p);
        disp('next step')
        disp(ie)
        disp(t_vals(end))
        tspan = [te(end):stepsize: p.endtime]; 
        c0 = [ye(end,1) ye(end,2)];
    
        t_vals_whole=vertcat(t_vals_whole,t_vals);
        c_vals_whole=vertcat(c_vals_whole,c_vals);
        r_vals=vertcat(r_vals,p.Rate*ones(length(t_vals),1));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        switch ie(end) 
            case 1 
                ie  %bug catching
                ie = [];
                "verified empty"%bug catching
                ie   %bug catching  
            case 2
                 p.Rate = p.Rate_0; 
            case 3 
                p.Rate=0;   
            otherwise
                
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isempty(ie)
            "verified empty"%bug catching
        end



    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dcdt = derivatives(t, c, p)
dcdt = [p.Rate/p.V1 - (p.k+p.k12)*c(1) + p.k21*p.V2/p.V1*c(2) 
        p.k12*p.V1/p.V2*c(1) - p.k21*c(2)                    ];          
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [value,isterminal,direction] = dose_events(t,c,p)
    %detect time to stop infusion 
    value(1,1) = t-p.endtime;
    isterminal(1,1) = 1;
    direction(1,1) = 1;

    %detect time to start infusion 
    value(2,1) = p.interval/2-mod(t,p.interval);
    isterminal(2,1) = 1;
    direction(2,1) = 1;
    
    %detect time to stop infusion 
    value(3,1) = mod(t,p.interval)-p.duration
    isterminal(3,1) = 1;
    direction(3,1) =1 ;

end