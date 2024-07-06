%synthetic data experiments with the Correlated Gaussian Wigner model.
% We compare MDGM,PDGGM and Grampa.
%   n_exp = identifier to the save the results of the
%   experiments.IMPORTANT: it overwrites if not changed
%   m_iter = maximum MD iterations
%   n = size of the graph
function []=experiments_wigner(n_exp,m_iter,n)
addpath(genpath('Helpers'));
addpath(genpath('Algos Graph Matching'));
num_run = 10; %number of Monte Carlo runs
vec_noise = 0:0.05:0.8; %grid of noise parameters
len_noise = length(vec_noise);
maxiter = m_iter;%maximum number of MD or GD iterations
%
corr_gr = zeros(len_noise, num_run);
corr_md = zeros(len_noise, num_run);
corr_pgd = zeros(len_noise, num_run);
%
run_gr = zeros(len_noise, num_run);
run_md = zeros(len_noise, num_run);
run_pgd = zeros(len_noise, num_run);
%
diag_gr = zeros(len_noise, num_run);
prop_gr = zeros(len_noise, num_run);
diag_md = zeros(len_noise, num_run);
prop_md = zeros(len_noise, num_run);
diag_pgd = zeros(len_noise, num_run);
prop_pgd = zeros(len_noise, num_run);
%% Iteration over independent samples 
for ind_run = 1:num_run
    fprintf('Iteration %i \n', ind_run);
        
        %% Iteration over noise levels 
        for ind_noise = 1:len_noise
            %random graph model
            sigma = vec_noise(ind_noise); disp(sigma);
             [A, B, P_rnd] = generate_wig_2(n,sigma,1);%generate wigner graphs
            % [A, B, P_rnd] = generate_wig_2(n,sigma,1);

            %% MD
            clock_md = tic;
            %initial=(1/n)*ones(n);
            [P,diag,prop] = matching_MD(A,B,maxiter,'dynamic_steps');
            fix_pt_ratio = sum(dot(P_rnd, P)) / n;
            corr_md(ind_noise, ind_run) = fix_pt_ratio;
            run_md(ind_noise, ind_run) = toc(clock_md);
            diag_md(ind_noise,ind_run) = diag;
            prop_md(ind_noise,ind_run) = prop;
                  
            %% Robust spectral method (aka Grampa)
            clock_grampa = tic;
            [P,diag,prop] = matching_robust_spectral(A, B, 0.2);
            fix_pt_ratio = sum(dot(P_rnd, P)) / n;
            corr_gr(ind_noise, ind_run) = fix_pt_ratio;
            run_gr(ind_noise, ind_run) = toc(clock_grampa);
            diag_gr(ind_noise,ind_run) = diag;
            prop_gr(ind_noise,ind_run) = prop;
                  
            
             %% PGD
            clock_pgd = tic;
            [P,diag,prop] = matching_PGD(A,B,maxiter,'heuristic');%,0.9);
            fix_pt_ratio = sum(dot(P_rnd, P)) / n;
            corr_pgd(ind_noise, ind_run) = fix_pt_ratio;
            run_pgd(ind_noise, ind_run) = toc(clock_pgd);
            diag_pgd(ind_noise,ind_run) = diag;
            prop_pgd(ind_noise,ind_run) = prop;
                  
            
        end
end

%compute the average overlap over the Monte Carlo runs
corr_md_mean = mean(corr_md, 2);
corr_gr_mean = mean(corr_gr, 2);
corr_pgd_mean=mean(corr_pgd,2);
%compute the average running time over the Monte Carlo runs
run_md_mean = mean(mean(run_md(1:end,:)));
run_gr_mean = mean(mean(run_gr(1:end,:)));
run_pgd_mean= mean(mean(run_pgd(1:end,:)));
%compute the average running time std the Monte Carlo runs
run_md_std = std(mean(run_md,2));
run_pgd_std = std(mean(run_pgd,2));
run_gr_std = std(mean(run_gr,2));

p1=0.95;
p2=0.05;
q1_pi_corr=quantile(corr_md,p1,2);
q1_pi2_corr=quantile(corr_gr,p1,2);
q1_pi3_corr=quantile(corr_pgd,p1,2);

q2_pi_corr=quantile(corr_md,p2,2);
q2_pi2_corr=quantile(corr_gr,p2,2);
q2_pi3_corr=quantile(corr_pgd,p2,2);

%disp([robust_run_mean,qm_PPM_corr_mean, full_qp_run_mean]);

    figure;hold on;
    hdata1=line(vec_noise, corr_md_mean);
    hdata2=line(vec_noise, corr_gr_mean);
    hdata3=line(vec_noise, corr_pgd_mean);
    curve11 = q1_pi_corr' ;
    curve12 = q2_pi_corr' ;
    hdata1_1 = [vec_noise, fliplr(vec_noise)];
    inBetween = [curve11, fliplr(curve12)];
    fill(hdata1_1, inBetween, [0.75 0 0],'FaceAlpha',0.2,'LineStyle',':','EdgeColor','None');
    curve21 = q1_pi2_corr' ;
    curve22 = q2_pi2_corr' ;
    hdata2_1 = [vec_noise, fliplr(vec_noise)];
    inBetween = [curve21, fliplr(curve22)];
    fill(hdata2_1, inBetween, [0.75 0 0.75],'FaceAlpha',0.2,'LineStyle',':','EdgeColor','None');
    curve31 = q1_pi3_corr';
    curve32 = q2_pi3_corr' ;
    hdata3_1 = [vec_noise, fliplr(vec_noise)];
    inBetween = [curve31, fliplr(curve32)];
    fill(hdata3_1, inBetween, [.75 .75 1],'FaceAlpha',0.2,'LineStyle',':','EdgeColor','None');
    set(hdata2, 'LineStyle','--','Marker', 'none', 'MarkerSize', 5,'MarkerEdgeColor', 'none','LineWidth', 2);
    set(hdata1, 'LineStyle','--','Color', [0.75 0 0],'Marker', 'none', 'MarkerSize', 5,'MarkerEdgeColor', 'r', 'MarkerFaceColor', [0.75 0 0],'LineWidth', 2);
    set(hdata3, 'LineStyle','--','Color', [.75 .75 1],'Marker', 'none', 'MarkerSize', 5,'MarkerEdgeColor', [.75 .75 1], 'MarkerFaceColor', [.75 .75 1],'LineWidth', 2);
    str = sprintf('Seedless GM comparison with n=%i, CGW',n);
    hTitle=title(str);set(hTitle, 'FontSize', 16, 'FontWeight' , 'bold')
    set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on','XColor', [.4 .4 .4], 'YColor', [.4 .4 .4],'LineWidth', 1,'Fontsize',14);
    %plot(vec_noise, m_pi_corr(j,:));hold on;plot(vec_noise, m_pi2_corr(j,:));hold on;plot(vec_noise, m_pi3_corr(j,:));hold on;plot(vec_noise, m_pi4_corr(j,:));
    hXLabel = xlabel('Noise $\sigma$','interpreter','latex');
    hYLabel = ylabel('Recovery fraction','interpreter','latex');
    set([hXLabel, hYLabel], 'FontSize', 20);
    
    hLegend = legend([hdata1, hdata2,hdata3],'EMDGM','Grampa','PGDGM');%,'PPMGM in.5','PPMGM in.6');%,'PPMGM in.5','PPMGM in.6');%,'Top-eigen','Top-eigen+PPM');
 
   set(hLegend, 'FontSize', 11);hold off

id_experiment=['/times_exp_' int2str(n_exp)];
destdirectory = strcat(pwd,'/Mat files');%folder where the data will be stored
id_savedvars = strcat(destdirectory,id_experiment);%identifier of the variables to be sote in the .mat file
%
save(id_savedvars,'run_md','run_pgd','run_gr','n','corr_md','corr_pgd','corr_gr','run_md_mean','run_gr_mean','run_pgd_mean','run_md_std','run_gr_std','run_pgd_std');%choose the variables to save
end 




