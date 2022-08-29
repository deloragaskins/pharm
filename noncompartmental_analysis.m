function noncompartmental_analysis()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters
% data genterating parameters
dg.k_elim= 0.5;
dg.ka=1
dg.c_0=15
dg.sampling_times=horzcat([15/60:15/60:1],[1.5:30/60:6]);
dg.drugamt=25
dg.noise_on=1
dg.mu=0
dg.sigma1=0.25
dg.sigma2=0.3
%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate data 
[t1,y1]= concentration_profile_IV(dg);
[t2,y2] = concentration_profile_PO(dg);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fit IV curve
curve_fit_params=IV_curve_fitter(t1,y1)
IVfit_t=[0:.1:8]
IVfit_y=exp(curve_fit_params.c_0+curve_fit_params.k_elim*IVfit_t)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plot results
fig_1=figure;
fig_1.Position = [200 200 1000 420];
ymax=max(max(y1),max(y2(:,2)))
fudge=1.2
x_last=dg.sampling_times(end)*1.5


%%%%%%%%%%%%%%
subplot(1,2,1)
plot(t1,y1,'or')
hold on 
semilogy(IVfit_t,IVfit_y,'-r')
%semilogy(t1,y1,'or')
xlim([0 x_last]);
ylim([-1 ymax*fudge]);

xlabel('\fontsize{13}Time [hours]')
ylabel('\fontsize{13} Concentration [mg/L]')
title( 'IV dose: Drug Concentration vs Time')
%%%%%%%%%%%%%%
subplot(1,2,2)
%plot(t2(:),y2(:,1),'xk')
plot(t2(:),y2(:,2),'dk')
xlim([0 x_last]);
ylim([-1 ymax*fudge]);
xlabel('\fontsize{13}Time [hours]')
ylabel('\fontsize{13}Concentration [mg/L]')
title( 'PO dose: Drug Concentration vs Time')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t,y] = concentration_profile_IV(params)
    y_0=params.c_0;
    t=params.sampling_times;
    y=y_0*exp(-params.k_elim*t);
    if params.noise_on==1
        rng(1,'twister');
        y=y+normrnd(params.mu,params.sigma1,size(y));
    end
end

function [t,y] = concentration_profile_PO(params)
    y_0 =[params.drugamt 0]; 
    options = odeset('RelTol',1e-12, 'AbsTol',[1e-12 1e-12]);
    tspan=params.sampling_times
    [t,y]= ode45(@derivatives, tspan, y_0, options, params);
    if params.noise_on==1
        rng(2,'twister');
        y=y+0.7*normrnd(params.mu,params.sigma2,size(y))
%         y(y<0)=0;
%         y
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fitted_params=IV_curve_fitter(t,y)
    %fitted_params.k_elim= 0.5;
    fitting_times=[0:0.2:t(end)];
    %reminder exp() is base e, base 10 is log10()
    
%fit attempts 1-3 failed, now i think its because I didnt keep good track
% of the exponential in the leading coefficient
%, but haven't gone back and checked. Should have 
% plotted y_transform vs t to catch this earlier 

% %%%%%%
% fit attempt #1
%     size(y_transform)
%     size(t)
%     X = [ones(size(t')) t'];
%     b = X \ y_transform;
%     fit_y=b(1)*exp(b(2)*fitting_times)
% %%%%%%
% fit attempt #2
%    y_transform=log(y)
%    times=t'
%    tbl = table(times,y_transform);
%    model1 = fitlm(tbl, 'y_transform ~ times')
% %%%%%%
% fit attempt #3
%     y_transform=log10(y);
%     p = polyfit(t,y,1);
%     fitline = 10.^(polyval(p,fitting_times));
%     p
% fit attempt #4 
    %intitial failure, installed curve fitting tool box  
    f1 = fit(t',y','exp1');
    fitline1=f1(fitting_times)
    y_transform=log(y)

    f2 = fit(t',y_transform','poly1')
    fitline2=exp(f2(fitting_times))
    fitted_params.c_0=f2.p2
    fitted_params.k_elim=f2.p1
    %%%%%%%%%%%%
    %plot
    fig_2=figure;
    fig_2.Position = [400 400 500 420];
    semilogy(t,exp(y_transform),'or')
    hold on 
    semilogy(fitting_times,fitline2,'ok')

    semilogy(t,y,'or')
    semilogy(fitting_times,fitline1,'ok')

    xlabel('\fontsize{13}Time [hours]')
    ylabel('\fontsize{13}Concentration [mg/L]')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydt = derivatives(t, y, params)
dydt = [-params.ka*y(1), 
         params.ka*y(1) - params.k_elim*y(2) ];
end

