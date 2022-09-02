%%%%%%%%%%%%%%%%%%%%%%%%
endtime=p.endtime
thresh_value=p.min_effective_conc
t_vals_1=[];
c_vals_1=[];
spacing_window=p.interval*3
unit_conv=24;
%%%%%%%%%%%%%%%%%%%
f = figure;
f.Position = [100 100 1050 400];
xlim_1=[0 endtime/unit_conv];
xtick_spacing=spacing_window/unit_conv;
ylim_1=[0 thresh_value*1.2]
plottime= linspace(0,endtime,100);

plot(t_vals_1/unit_conv,c_vals_1(:,2),'-.b','DisplayName','conc of interest')
hold on
plot(plottime,thresh_value*ones(size(plottime)),'-c','DisplayName','Threshold Value')

xlim(xlim_1);
xticks(xlim_1(1):xtick_spacing:xlim_1(2))
ylim(y_lim_1)
title( 'Title')
legend('Location','southeast')
xlabel('\fontsize{13}Time [units]')
ylabel('\fontsize{13}Concentration [units]')
%set(gca,"FontSize",10)

saveas(gcf,'filename.png')