function noncompartmental_analysis()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters
% data Genterating parameters
dg.k_elim= 0.5;
dg.ka=1
dg.c_0=15
dg.sampling_times=horzcat([15/60:15/60:1],[1.5:30/60:6]);
dg.drugamt=25
dg.noise_on=1
dg.mu=0
dg.sigma=.15

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate data 
[t1,y1]= concentration_profile_IV(dg);
[t2,y2] = concentration_profile_PO(dg);
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
%semilogy(t1,y1,'or')
xlim([0 x_last]);
ylim([0 ymax*fudge]);

xlabel('\fontsize{13}Time [hours]')
ylabel('\fontsize{13} Concentration [mg/L]')
title( 'IV dose: Drug Concentration vs Time')
%%%%%%%%%%%%%%
subplot(1,2,2)
%plot(t2(:),y2(:,1),'xk')
plot(t2(:),y2(:,2),'dk')
xlim([0 x_last]);
ylim([0 ymax*fudge]);
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
        y=y+normrnd(params.mu,params.sigma,size(y));
    end
end

function [t,y] = concentration_profile_PO(params)
    y_0 =[params.drugamt 0]; 
    options = odeset('RelTol',1e-12, 'AbsTol',[1e-12 1e-12]);
    tspan=params.sampling_times;
    [t,y]= ode45(@derivatives, tspan, y_0, options, params);
    if params.noise_on==1
        rng(2,'twister');
        y=y+normrnd(params.mu,params.sigma,size(y))
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydt = derivatives(t, y, params)
dydt = [-params.ka*y(1), 
         params.ka*y(1) - params.k_elim*y(2) ];
end

