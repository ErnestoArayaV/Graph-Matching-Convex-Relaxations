function plot_autonomous_results(curve_md,curve_pgd,curve_sp,str_ylabel,curve_opt)
if (~exist('curve_opt', 'var'))
        logical = false;
else 
    logical = true;
end
figure;hold on;
dates =["03/31","04/07","04/14","04/21","04/28","05/05","05/12","05/19","05/26"];
cat = categorical(dates);
hdata1=line(cat, curve_md);
hdata2=line(cat, curve_sp);
hdata3=line(cat, curve_pgd);

set(hdata2, 'LineStyle','--','Marker', 'x', 'MarkerSize', 5,'MarkerEdgeColor', 'b','LineWidth', 2);
set(hdata1, 'LineStyle','--','Color', [0.75 0 0],'Marker', 'o', 'MarkerSize', 5,'MarkerEdgeColor', 'r', 'MarkerFaceColor', [0.75 0 0],'LineWidth', 2);
set(hdata3, 'LineStyle','--','Color', [.75 .75 1],'Marker', '+', 'MarkerSize', 5,'MarkerEdgeColor', [.75 .75 1], 'MarkerFaceColor', [.75 .75 1],'LineWidth', 2);
%hTitle=title(str);set(hTitle, 'FontSize', 16, 'FontWeight' , 'bold')
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on','XColor', [.4 .4 .4], 'YColor', [.4 .4 .4],'LineWidth', 1,'Fontsize',14);
%plot(vec_noise, m_pi_corr(j,:));hold on;plot(vec_noise, m_pi2_corr(j,:));hold on;plot(vec_noise, m_pi3_corr(j,:));hold on;plot(vec_noise, m_pi4_corr(j,:));
hXLabel = xlabel('Time','interpreter','latex');
%hYLabel = ylabel('Recovery fraction','interpreter','latex');
hYLabel = ylabel(str_ylabel,'interpreter','latex');
if logical
    hdata4=line(cat, curve_opt);
    set(hdata4, 'LineStyle','--','Color', [0.9290 0.6940 0.1250],'Marker', '^', 'MarkerSize', 5,'MarkerEdgeColor', [0.9290 0.6940 0.1250], 'MarkerFaceColor', [0.9290 0.6940 0.1250],'LineWidth', 2);
    hLegend = legend([hdata1, hdata2,hdata3,hdata4],'EMDGM','Grampa','PGDGM','Truth');
else 
    hLegend = legend([hdata1, hdata2,hdata3],'EMDGM','Grampa','PGDGM');
end 
set([hXLabel, hYLabel], 'FontSize', 20);

%hLegend = legend([hdata1, hdata2,hdata3],'EMDGM','Grampa','PGDGM');
set(hLegend, 'FontSize', 11);hold off


