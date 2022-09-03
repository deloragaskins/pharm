  
function [t_vals_whole, c_vals_whole]=run_dosing_course(derivativefunction,t_sim_end,dosing_interval,initial_dose,dose,p)   
    %set up integration
    t_vals_whole=[];
    c_vals_whole=[];
    tspan = [0 dosing_interval];     % max. time domain (h)
    c0 =[initial_dose,0, 0];
    options = odeset('RelTol',1e-12, 'AbsTol',[1e-12 1e-12 1e-12]);
    %run stepwise integration
    while tspan(1)<t_sim_end
        [t_vals,c_vals]= ode45(derivativefunction, tspan, c0, options, p);
        tspan = [t_vals(end) t_vals(end)+dosing_interval]; 
        c0 = [c_vals(end,1)+dose, c_vals(end,2), c_vals(end,3)];
        t_vals_whole=vertcat(t_vals_whole,t_vals);
        c_vals_whole=vertcat(c_vals_whole,c_vals);
    end
end
