function [dose_interval, dose_amount,max_dose] = dosing_chooser(derivativefunction, min_effective_conc, max_tolerated_conc,deriv_params)    
    %required parameters 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %for @derivatives
    %p.ka
    %p.k_elim
    %p.k12
    %p.k21
    %p.V1
    %p.V2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %for dosing_chooser
    %p.min_effective_conc
    %p.max_tolerated_conc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t_end_1=72;
    t_end_2=t_end_1*15;
    tspan1 = [0 t_end_1];     % max. time domain (h)
    f_0 = figure;
    unit_conv=24;
    %xlim_1=[0 t_end_1/unit_conv];
    xlim_2=[0 t_end_2/unit_conv];
    plottime= linspace(0, t_end_1,100);  
    plottime2= linspace(0,t_end_2,100);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %find max dose
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     doses=[0 10 100 1000 10000];
%     doses=[1 2 4 6 8 10]*1000;
%     doses=2000+ [0 1 2 3 4 5 6 7 8 9 10 ]*100;
%     doses=2100+[1 2 3 4 5 6 7 8 9 10]*10;
    doses=2150;
    %%%%%%%%%%%%%%%%%%%
    %sweeping for this parameteclosr 
    max_dose=2150;
    %%%%%%%%%%%%%%%%
    for initial_dose = doses
        c0 =[initial_dose,0, 0];
        options = odeset('RelTol',1e-12, 'AbsTol',[1e-12 1e-12 1e-12]);
        [t_vals,c_vals]= ode45(derivativefunction, tspan1, c0, options, deriv_params);
        
        plot (t_vals,c_vals(:,2),'-k')
        hold on
        plot (t_vals,c_vals(:,3),'-b')
        plot (t_vals,c_vals(:,1),'-r')

    end
    plot(plottime,(max_tolerated_conc)*ones(size(plottime)),'-c','DisplayName','threshold span')
    plot(plottime,(min_effective_conc)*ones(size(plottime)),'-c','DisplayName','threshold span')
    ylim([0 1.1*max_tolerated_conc])
    xlabel('\fontsize{13}Time [hours]')
    ylabel('\fontsize{13}Concentration [mg/L]')
    title( 'Single Dose to find max dose')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %find combination of dose amount/dose interval which keeps concentration
    %in between max tolterated concentration and min effective dose
    f_1 = figure;

    %%%%%%%%%%%%%%%%%%%%%
    %sweeping for these params
    dose_interval=24/2;
    dose_amount=max_dose/5;
    %%%%%%%%%%%%%%%%%%%%%%%%

    [t_vals_whole, c_vals_whole]=run_dosing_course(derivativefunction,t_end_2,dose_interval,dose_amount,dose_amount,deriv_params); 

    plot (t_vals_whole/unit_conv,c_vals_whole(:,2),'-k')
    hold on
    plot (t_vals_whole/unit_conv,c_vals_whole(:,3),'-b')

    plot(plottime2,(max_tolerated_conc)*ones(size(plottime2)),'-c','DisplayName','MTC')
    plot(plottime2,(min_effective_conc)*ones(size(plottime2)),'-c','DisplayName','MEC')
    
    ylim([0 1.1*max_tolerated_conc])
    xlim(xlim_2)
    xlabel('\fontsize{13}Time [days]')
    ylabel('\fontsize{13}Concentration [mg/L]')
    title( 'Dosing untill SS is achieved')

    
end 