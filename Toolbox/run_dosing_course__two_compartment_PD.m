function [t_vals_whole, y_vals_whole]=run_dosing_course__two_compartment_PD(derivativefunction,t_sim_end,dosing_interval,initial_dose,dose,initial_n,deriv_params)   
    %set up integration
    t_vals_whole=[];
    y_vals_whole=[];
    tspan = [0 dosing_interval];     % max. time domain (h)
    y0 =[initial_dose,0, 0,initial_n];
    options = odeset('RelTol',1e-12, 'AbsTol',[1e-12 1e-12 1e-12 1e-12]);
    %run stepwise integration
    while tspan(1)<t_sim_end
        [t_vals,y_vals]= ode15s(derivativefunction, tspan, y0, options, deriv_params);
        tspan = [t_vals(end) t_vals(end)+dosing_interval]; 
        y0 = [y_vals(end,1)+dose, y_vals(end,2), y_vals(end,3), y_vals(end,4)];
        t_vals_whole=vertcat(t_vals_whole,t_vals);
        y_vals_whole=vertcat(y_vals_whole,y_vals);
       
    end
end
