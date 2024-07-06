%function to simulate and save experiments for autonomous systems graphs.
%Based on code by Jiaming Xu for "Spectral Graph Matching and Regularized Quadratic Relaxations I: The Gaussian Model"
function [] = auto_sys_dynamics(n_exp,n)
addpath(genpath('Helpers'));
addpath(genpath('Algos Graph Matching'));
path = './Auto sys/';
load1=load(strcat(path,'/mat_files/auto_sys_data.mat'));
load2=load(strcat(path,'/mat_files/grad_norms_n1000_it10000.mat'));
auto_sys_mat=load1.auto_sys_mat;
norms_inf = load2.norms_inf;
norms_sq = load2.norms_sq;
%n = 1000; %number of samples
deg = sum(auto_sys_mat{1});
[~, I] = maxk(deg, n);
A = auto_sys_mat{1}(I, I); %creatte first graph
%test standardization
%A = standardization_empirical(A); %uncomment to standardize% see our
%paper for details

max_iter = (log(n)^3);%maximum number of iterations for MD and PGD
%
%metrics Grampa
corr_sp = zeros(9, 1);
val_sp = zeros(9, 1);
%
%metrics PGD
corr_pgd4 = zeros(9, 1);
val_pgd4 = zeros(9, 1);
%rates_pgd4 = zeros(9,max_iter);
%
%metrics MD
corr_md2 = zeros(9, 1);
val_md2 = zeros(9, 1);
%rates_md2 = zeros(9,max_iter);
% 
val_truth = zeros(9, 1);
val_cnv = zeros(9,1);
%
%runing times
rtime_sp = zeros(9,1);
rtime_pgd4 = zeros(9,1);
rtime_md = zeros(9,1);
%
for i = 1:9
    disp(i);
    B = auto_sys_mat{i}(I, I);
    %test standardization
    %B = standardization_empirical(B);
    val_truth(i) = sum(dot(A, B));
    val_cnv(i) = convex_obj(A, B, speye(n));
    
    P_rnd = sparse(eye(n));
    P_rnd = P_rnd(:, randperm(n));
    B = P_rnd * B * P_rnd';

    %% Robust spectral (grampa)
    clock_grampa = tic;
  %  [P_sp,X_grampa] = matching_robust_spectral(full(A), full(B), 1);
    [P_sp,~] = matching_robust_spectral(full(A), full(B), 1);
    P_sp=sparse(P_sp);
    corr_sp(i) = sum(dot(P_rnd, P_sp))/n;
    val_sp(i) = sum(dot(P_sp * A * P_sp', B));
    rtime_sp(i) = toc(clock_grampa);
    %norm_X_grampa = X_grampa/sum(sum(abs(X_grampa)));
    %obj_grampa = convex_obj(A,B,norm_X_grampa);
%         
     %% MD2
    clock_md2 = tic;
    P_md = matching_MD(sparse(A),sparse(B),max_iter,'fixed_steps',norms_inf(i));
    corr_md2(i) = sum(dot(P_rnd, P_md))/n;
    val_md2(i) = sum(dot(P_md * A * P_md', B));
    rtime_md(i) = toc(clock_md2);

    %% PGD 4
    clock_pgd4 = tic;
    P_pgd = sparse(matching_PGD(sparse(A),sparse(B),max_iter,'fixed_steps',norms_sq(i)));
    corr_pgd4(i) = sum(dot(P_rnd, P_pgd))/n;
    val_pgd4(i) = sum(dot(P_pgd * A * P_pgd', B));
    rtime_pgd4(i) = toc(clock_pgd4);
        
    clear B P_rnd P_sp P_md;
end
id_experiment=[ '/comparison_exp_' int2str(n_exp)];
destdirectory = strcat(pwd,'/Mat files');%folder where the data will be stored
id_savedvars = strcat(destdirectory,id_experiment);%identifier of the variables to be sote in the .mat file
%save('.\mat_files\auto_sys_dynamics.mat');
save(id_savedvars, 'corr_sp', 'corr_md2', 'corr_pgd4','val_sp', 'val_md2','val_pgd4','val_truth');