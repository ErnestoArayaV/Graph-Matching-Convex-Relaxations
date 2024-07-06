function plot_nice_fb(s,data,p1,p2)
%compute mean over Monte Carlos
data_m = mean(data,3);
%compute quantiles
q1_data1=quantile(data(1,:,:),p1,3);
q1_data2=quantile(data(2,:,:),p1,3);
q1_data3=quantile(data(3,:,:),p1,3);
%
q2_data1=quantile(data(1,:,:),p2,3);
q2_data2=quantile(data(2,:,:),p2,3);
q2_data3=quantile(data(3,:,:),p2,3);
%plot
figure;hold on;
% plot(s,correct_m(2,:),'+--','Color',[0.75 0 0]);%plot md
% plot(s,correct_m(1,:),'o--','Color','blue');%plot grampa
% plot(s,correct_m(3,:),'^--','Color',[.75 .75 1]);%plot pgd
hdata1=line(s, data_m(2,:));%grampa
hdata2=line(s, data_m(1,:));%md
hdata3=line(s, data_m(3,:));%pgd
curve11 = q1_data1 ;
curve12 = q2_data1 ;
hdata1_1 = [s, fliplr(s)];
inBetween = [curve11, fliplr(curve12)];
fill(hdata1_1, inBetween, [0.75 0 0.75],'FaceAlpha',0.2,'LineStyle',':','EdgeColor','None');
curve21 = q1_data2 ;
curve22 = q2_data2 ;
hdata2_1 = [s, fliplr(s)];
inBetween = [curve21, fliplr(curve22)];
fill(hdata2_1, inBetween,[0.75 0 0] ,'FaceAlpha',0.2,'LineStyle',':','EdgeColor','None');
curve31 = q1_data3;
curve32 = q2_data3 ;
hdata3_1 = [s, fliplr(s)];
inBetween = [curve31, fliplr(curve32)];
fill(hdata3_1, inBetween, [.75 .75 1],'FaceAlpha',0.2,'LineStyle',':','EdgeColor','None');
set(hdata2, 'LineStyle','--','Marker', '+', 'MarkerSize', 5,'MarkerEdgeColor', 'blue','LineWidth', 2);
set(hdata1, 'LineStyle','--','Color', [0.75 0 0],'Marker', 'o', 'MarkerSize', 5,'MarkerEdgeColor', 'r', 'MarkerFaceColor', [0.75 0 0],'LineWidth', 2);
set(hdata3, 'LineStyle','--','Color', [.75 .75 1],'Marker', '^', 'MarkerSize', 5,'MarkerEdgeColor', [.75 .75 1], 'MarkerFaceColor', [.75 .75 1],'LineWidth', 2);
%str = sprintf('Seedless GM comparison with n=%i, CER',n);
%hTitle=title(str);set(hTitle, 'FontSize', 16, 'FontWeight' , 'bold')
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on','XColor', [.4 .4 .4], 'YColor', [.4 .4 .4],'LineWidth', 1,'Fontsize',14);
%plot(vec_noise, m_pi_corr(j,:));hold on;plot(vec_noise, m_pi2_corr(j,:));hold on;plot(vec_noise, m_pi3_corr(j,:));hold on;plot(vec_noise, m_pi4_corr(j,:));
hXLabel = xlabel('$s$','interpreter','latex');
hYLabel = ylabel('Recovery fraction','interpreter','latex');
set([hXLabel, hYLabel], 'FontSize', 20);

hLegend = legend([hdata1, hdata2,hdata3],'EMDGM','Grampa','PGDGM');%,'PPMGM in.5','PPMGM in.6');%,'PPMGM in.5','PPMGM in.6');%,'Top-eigen','Top-eigen+PPM');

set(hLegend, 'FontSize', 11);hold off

