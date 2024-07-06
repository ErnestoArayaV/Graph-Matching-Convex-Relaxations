% Compare seedless GM methods. Grampa, Umeyama, full QP and Mirror descent.
clear all;
addpath(genpath('Helpers'));
addpath(genpath('Algos Graph Matching'));

%% initialization
p = 0.5;
num_run = 15;
vec_dim = 300;%500:10:500;
len_dim = length(vec_dim);
vec_noise = 0:0.05:0.8; %grid of noise parameters
len_noise = length(vec_noise);

full_qp_corr = ones(len_dim, len_noise, num_run);% to record overlap with gt in full QP
full_qp_run = zeros(len_dim, len_noise, num_run);% to record running time in full QP

robust_corr = ones(len_dim, len_noise, num_run);% to record overlap with gt in Grampa
robust_run = zeros(len_dim, len_noise, num_run);% to record running time in Grampa

MD_simplex_corr=ones(len_dim, len_noise, num_run);% to record overlap with gt in simplex MD
MD_simplex_run=zeros(len_dim, len_noise, num_run);% to record running time in in simplex MD

umeyama_corr=ones(len_dim, len_noise, num_run);% to record overlap with gt in umeyama
umeyama_run=zeros(len_dim, len_noise, num_run);% to record running time in in umeyama

%% Iteration over independent samples 
for ind_run = 1:num_run
    fprintf('Iteration %i \n', ind_run);

    %% Iteration over dimensions 
    for ind_dim = 1:len_dim
        n = vec_dim(ind_dim);
        fprintf('Matrix dimension %i \n', n);
        maxiter = (log(n)^3);
        
        %% Iteration over noise levels 
        for ind_noise = 3:len_noise
            %correlated random graph model
            sigma = vec_noise(ind_noise); disp(sigma);
            [A, B, A0, B0, P_rnd] = generate_er(n, p, sigma);
            % [A, B, P_rnd] = generate_wig_2(n,sigma,1);
            %% Full QP 
            tic;
            P = matching_full_qp(A, B);
            fix_pt_ratio = sum(dot(P_rnd, P)) / n;
            full_qp_corr(ind_dim, ind_noise, ind_run) = fix_pt_ratio;
            full_qp_run(ind_dim, ind_noise, ind_run) = toc;
            %% MD
            tic;
            %initial=ones(n)/(n^2);
            [P,obj_MD] = matching_MD(A,B,maxiter,'dynamic_steps');
            %[P,~,~,~,~,~,~,~,~,~,~,~]=matching_simplex_mirrordescent(A,B,maxiter,initial/n,P_rnd);
            fix_pt_ratio = sum(dot(P_rnd, P)) / n;
            MD_simplex_corr(ind_dim, ind_noise, ind_run) = fix_pt_ratio;
            MD_simplex_run(ind_dim, ind_noise, ind_run) = toc;
                  
            %% Robust spectral method (Grampa)
            tic;
            P = matching_robust_spectral(A, B, 0.2);
            fix_pt_ratio = sum(dot(P_rnd, P)) / n;
            robust_corr(ind_dim, ind_noise, ind_run) = fix_pt_ratio;
            robust_run(ind_dim, ind_noise, ind_run) = toc;   

            %% Umeyama
            tic;
             P = matching_umeyama(A, B);
            fix_pt_ratio = sum(dot(P_rnd, P)) / n;
            umeyama_corr(ind_dim, ind_noise, ind_run) = fix_pt_ratio;
            umeyama_run(ind_dim, ind_noise, ind_run) = toc;  
            
        end
    end
end
%compute the average overlap over the Montecarlo runs
full_qp_corr_mean = mean(full_qp_corr, 3);
robust_corr_mean = mean(robust_corr, 3);
MD_simplex_corr_mean=mean(MD_simplex_corr,3);
umeyama_corr_mean = mean(umeyama_corr,3); 

%compute the average running time over the Montecarlo runs
full_qp_run_mean = mean(mean(mean(full_qp_run(1, 3:end, :))));
robust_run_mean = mean(mean(mean(robust_run(1, 3:end, :))));
MD_simplex_run_mean=mean(mean(mean(MD_simplex_run(1, 3:end, :))));
umeyama_run_mean=mean(mean(mean(umeyama_run(1, 3:end, :))));

%compute the quantiles for the 'confidence' regions
p1=0.95;
p2=0.05;
q1_pi3_corr=quantile(full_qp_corr,p1,3);
q1_pi_corr=quantile(MD_simplex_corr,p1,3);
q1_pi2_corr=quantile(robust_corr,p1,3);
q1_pi4_corr=quantile(umeyama_corr,p1,3);

q2_pi3_corr=quantile(full_qp_corr,p2,3);
q2_pi_corr=quantile(MD_simplex_corr,p2,3);
q2_pi2_corr=quantile(robust_corr,p2,3);
q2_pi4_corr=quantile(umeyama_corr,p2,3);

for j=1:len_dim
    figure;hold on;
    hdata3=line(vec_noise, full_qp_corr_mean(j,:));
    hdata1=line(vec_noise, MD_simplex_corr_mean(j,:));
    hdata2=line(vec_noise, robust_corr_mean(j,:));
    hdata4 = line(vec_noise, umeyama_corr_mean(j,:));
    
    curve11 = q1_pi_corr (j,:);
    curve12 = q2_pi_corr (j,:);
    hdata1_1 = [vec_noise, fliplr(vec_noise)];
    inBetween = [curve11, fliplr(curve12)];
    fill(hdata1_1, inBetween, [0.75 0 0],'FaceAlpha',0.2,'LineStyle',':','EdgeColor','None');
    curve21 = q1_pi2_corr (j,:);
    curve22 = q2_pi2_corr (j,:);
    hdata2_1 = [vec_noise, fliplr(vec_noise)];
    inBetween = [curve21, fliplr(curve22)];
    fill(hdata2_1, inBetween, [0.75 0 0.75],'FaceAlpha',0.2,'LineStyle',':','EdgeColor','None');
    curve31 = q1_pi3_corr (j,:);
    curve32 = q2_pi3_corr (j,:);
    hdata3_1 = [vec_noise, fliplr(vec_noise)];
    inBetween = [curve31, fliplr(curve32)];
    fill(hdata3_1, inBetween, [0.4660 0.6740 0.1880],'FaceAlpha',0.2,'LineStyle',':','EdgeColor','None');
    curve41 = q1_pi4_corr (j,:);
    curve42 = q2_pi4_corr (j,:);
    hdata4_1 = [vec_noise, fliplr(vec_noise)];
    inBetween = [curve41, fliplr(curve42)];
    fill(hdata4_1, inBetween, [0.9290 0.6940 0.1250],'FaceAlpha',0.2,'LineStyle',':','EdgeColor','None');
    
    set(hdata2, 'LineStyle','--','Marker', 'none', 'MarkerSize', 5,'MarkerEdgeColor', 'none','LineWidth', 2);
    set(hdata1, 'LineStyle','--','Color', [0.75 0 0],'Marker', 'none', 'MarkerSize', 5,'MarkerEdgeColor', 'r', 'MarkerFaceColor', [0.75 0 0],'LineWidth', 2);
    set(hdata3, 'LineStyle','--','Color', [0.4660 0.6740 0.1880],'Marker', 'none', 'MarkerSize', 5,'MarkerEdgeColor', [0.4660 0.6740 0.1880], 'MarkerFaceColor', [0.4660 0.6740 0.1880],'LineWidth', 2);
    set(hdata4, 'LineStyle','--','Color', [0.9290 0.6940 0.1250],'Marker', 'none', 'MarkerSize', 5,'MarkerEdgeColor', [0.9290 0.6940 0.1250], 'MarkerFaceColor', [0.9290 0.6940 0.1250],'LineWidth', 2);

    str = sprintf('Seedless GM comparison with n=%i, CER',vec_dim(j));
    hTitle=title(str);set(hTitle, 'FontSize', 16, 'FontWeight' , 'bold')
    set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on','XColor', [.4 .4 .4], 'YColor', [.4 .4 .4],'LineWidth', 1,'Fontsize',14);
    hXLabel = xlabel('Noise $\sigma$','interpreter','latex');
    hYLabel = ylabel('Recovery fraction','interpreter','latex');
    set([hXLabel, hYLabel], 'FontSize', 20);
    
    hLegend = legend([hdata1, hdata2,hdata3,hdata4],'EMDGM','Grampa','QPADMM','Umeyama');%,'PPMGM in.5','PPMGM in.6');%,'PPMGM in.5','PPMGM in.6');%,'Top-eigen','Top-eigen+PPM');
 
    set(hLegend, 'FontSize', 11);hold off
end

%clear -regexp _corr$ _run$;
id_experiment='/comparison_exp_3';% IMPORTANT: this has to be changed if we do not want to overwrite the files. 
destdirectory = strcat(pwd,'/Mat files');%folder where the data will be stored
id_savedvars = strcat(destdirectory,id_experiment);%identifier of the variables to be sote in the .mat file
save(id_savedvars,'robust_corr');
save(id_savedvars,'MD_simplex_corr',"-append");
save(id_savedvars,'full_qp_corr',"-append");
save(id_savedvars,'robust_run',"-append");
save(id_savedvars,'full_qp_run',"-append");
save(id_savedvars,'MD_simplex_run',"-append");


