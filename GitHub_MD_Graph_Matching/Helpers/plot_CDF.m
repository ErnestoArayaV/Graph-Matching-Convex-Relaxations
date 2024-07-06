%function the plotting the CDF for the SHREC 16 computer vision dataset

function plot_CDF(thr,curve1,curve2,opt)
if exist('opt','var')
    curve3 = opt.curve;
    curve3 = curve3/100;
    curve3_m = mean(curve3,1);
end
curve1=curve1/100;
curve2 = curve2/100;
%compute the average
curve1_m = mean(curve1, 1);
curve2_m = mean(curve2, 1);
%
p1=0.95;
p2=0.05;
%
q1_data1=quantile(curve1,p1,1);
q1_data2=quantile(curve2,p1,1);
q1_data3=quantile(curve3,p1,1);
%
q2_data1=quantile(curve1,p2,1);
q2_data2=quantile(curve2,p2,1);
q2_data3=quantile(curve3,p2,1);
%
figure;hold on;
%thr=thr';

hdata1=line(thr, curve1_m);
hdata2=line(thr, curve2_m);
if exist('opt','var')
    if opt.plot_std == 1
        
    curve11 = q1_data1 ;
    curve12 = q2_data1 ;
    hdata1_1 = [thr, fliplr(thr)];
    inBetween = [curve11, fliplr(curve12)];
    fill(hdata1_1, inBetween, [0.75 0 0],'FaceAlpha',0.2,'LineStyle',':','EdgeColor','None');
    curve21 = q1_data2 ;
    curve22 = q2_data2 ;
    hdata2_1 = [thr, fliplr(thr)];
    inBetween = [curve21, fliplr(curve22)];
    fill(hdata2_1, inBetween, [0.75 0 0.75],'FaceAlpha',0.2,'LineStyle',':','EdgeColor','None');
    end
end

set(hdata2, 'LineStyle','-','Marker', 'none', 'MarkerSize', 5,'MarkerEdgeColor', 'none','LineWidth', 2);
set(hdata1, 'LineStyle','-','Color', [0.75 0 0],'Marker', 'none', 'MarkerSize', 5,'MarkerEdgeColor', 'r', 'MarkerFaceColor', [0.75 0 0],'LineWidth', 2);
if exist('opt','var')
    hdata3=line(thr, curve3_m);
    curve31 = q1_data3;
    curve32 = q2_data3 ;
    %hdata3_1 = [thr, fliplr(thr)];
    set(hdata3, 'LineStyle','-','Color', [.75 .75 1],'Marker', 'none', 'MarkerSize', 5,'MarkerEdgeColor', [.75 .75 1], 'MarkerFaceColor', [.75 .75 1],'LineWidth', 2);
    if opt.plot_std == 1
        
        hdata3_1 = [thr, fliplr(thr)];
        inBetween = [curve31, fliplr(curve32)];
        %set(hdata3, 'LineStyle','--','Color', [.75 .75 1],'Marker', 'none', 'MarkerSize', 5,'MarkerEdgeColor', [.75 .75 1], 'MarkerFaceColor', [.75 .75 1],'LineWidth', 2);
        fill(hdata3_1, inBetween, [.75 .75 1],'FaceAlpha',0.2,'LineStyle',':','EdgeColor','None');
    end
    hLegend = legend([hdata1, hdata2,hdata3],'EMDGM','Grampa','PGDGM');
else
    hLegend = legend([hdata1, hdata2],'EMDGM','Grampa');
end

%str = sprintf('Seedless GM comparison with n=%i, CGW',n);
%hTitle=title(str);set(hTitle, 'FontSize', 16, 'FontWeight' , 'bold')
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on','XColor', [.4 .4 .4], 'YColor', [.4 .4 .4],'LineWidth', 1,'Fontsize',14);
%plot(vec_noise, m_pi_corr(j,:));hold on;plot(vec_noise, m_pi2_corr(j,:));hold on;plot(vec_noise, m_pi3_corr(j,:));hold on;plot(vec_noise, m_pi4_corr(j,:));
hXLabel = xlabel('$\epsilon$','interpreter','latex');
hYLabel = ylabel('CDF','interpreter','latex');
set([hXLabel, hYLabel], 'FontSize', 20);

%hLegend = legend([hdata1, hdata2,hdata3],'EMDGM','Grampa','PGDGM');%,'PPMGM in.5','PPMGM in.6');%,'PPMGM in.5','PPMGM in.6');%,'Top-eigen','Top-eigen+PPM');

set(hLegend, 'FontSize', 11);hold off