% comparison emdgm, pgdgm and grampa algorithms. Based on a code by  
% applied to the SHREC16 dataset

%% PREPARE
clear all;
path_kids = './TOPKIDS';          % path to the complete TOPKIDS data set
track = 'low resolution/';  % low or high resolution
addpath(genpath('Helpers'));
addpath(genpath('Algos Graph Matching'));
%% CALCULATE CURVES

thresholds = 0:0.01:1;

nb_examples = 90;
curve1 = zeros(nb_examples,length(thresholds));
curve2 = zeros(nb_examples,length(thresholds));
for k=1:min(nb_examples,90)
    k
    i=ceil(k/9);
    j=k-(i-1)*9;
    if j<i
        j=j+15;
    else
        j=j+16;
    end
    i=i+15;
    %% load data
    
    M = load_off(strcat(path_kids,track, 'kid', num2str(i,'%02d'), '.off'));
    N = load_off(strcat(path_kids,track, 'kid', num2str(j,'%02d'), '.off'));
    
    V1=M.VERT;                  % 3-d coordinates of vertices
    F1=M.TRIV;                  % face for triangulation
    V2=N.VERT;
    F2=N.TRIV;
    
    adj1 = triangulation2adjacency(M.TRIV);     % adj after triangulation
    adj2 = triangulation2adjacency(N.TRIV);
    dist=geodesic_distance(N.TRIV,N.VERT); %Added by JX
    dist=sparse(dist);
    n1=size(adj1,1);
    n2=size(adj2,1);
    
    EYE1=sparse(1:n1,1:n1,1,n1,n1);
    EYE2=sparse(1:n2,1:n2,1,n2,n2);
    W21=double((double((adj1*adj1)>0)-adj1-EYE1)>0);
    W22=double((double((adj2*adj2)>0)-adj2-EYE2)>0);
%% Read ground truth
    gt_M_null = read_correspondence(strcat(path_kids, track, 'kid', num2str(i,'%02d'), '_ref.txt'));
    gt_N_null = read_correspondence(strcat(path_kids, track, 'kid', num2str(j,'%02d'), '_ref.txt'));
    gt = merge_ground_truth(gt_M_null, gt_N_null);
    P_rnd=zeros(n2,n1);
    for ind=1:length(gt(:,1))
        P_rnd(gt(ind,2),gt(ind,1))=1;
    end
    P_rnd=sparse(P_rnd);

%% MD algorithm
    largest_grad_norm =0;
    lar_gradnorm =0;
    %curve1 = zeros(1,length(thresholds));
    iter_max=100;
    %r_old = 1;
    %r_old=corr_sp;
    adj1 = adj1/full(sqrt(sum(abs(adj1(:)).^2)));
    adj2 = adj2/full(sqrt(sum(abs(adj2(:)).^2)));
    %S=P_sp';
    adj2_sq = sparse(adj2^2);
    adj1_sq = sparse(adj1^2);
    X = ones(size(adj2,1),size(adj1,1))/(n1*n2);
    X2 =X;
    learning_rate = 2.91e+10;%fixed learning rate. Alternatively, it can be defined to be decreasing inside the loop
    obj_best_1=1000000;
    obj_best_2=1000000;
    %MD loop and PGD loops
    for iter_count=1:1:iter_max
        %print iterations 
        if mod(iter_count,20) == 0
            iter_count
            %learning_rate = learning_rate/(iter_count/5);%define here the
            %decreasing learning rate
        end
        %MD explicit multiplicative updates
        grad = 2*(adj2_sq*X+X*adj1_sq)-4*adj2*X*adj1;
        grad_norm = max(max(abs(grad)));
        gradnorm = norm(grad,'fro');
        lar_gradnorm = max(gradnorm,lar_gradnorm);
        largest_grad_norm = max([grad_norm,largest_grad_norm]);
        if grad_norm<0.0000001
            learning_rate = 2*sqrt(2)/(largest_grad_norm*sqrt(iter_count+1));
        else 
            learning_rate = 2*sqrt(2)/(grad_norm*(iter_count+1)^.5);
        end 
         actual_gamma=(0.7*gradnorm^2)/norm(adj2*grad-grad*adj1,'fro')^2;
         X = (X.*exp(-learning_rate*grad))/sum(sum(X.*exp(-learning_rate*grad))); %MD update
         X2 = X2-actual_gamma*grad;%PGD update
         x2=reshape(X2, [], 1);
         y=projsplx(x2);
         X2=reshape(y,[n2,n1]);
        %Rounding step. We do it here to have a convergence criterion
        actual_obj_1 = norm(adj2*X-X*adj1,'fro');
        if actual_obj_1<obj_best_1
            S_best=S1;
            obj_best_1 =actual_obj_1;
        end
        actual_obj_2 = norm(adj2*X2-X2*adj1,'fro');
        if actual_obj_2<obj_best_2
            S_best2=S2;
            obj_best_2 =actual_obj_2;
        end
        
    end
    S=S_best;
    S2=S_best2;
  
    if n1<=n2
        final_idx1 = [1:n1]';
        idx2=1:n2;
        final_idx2 = S'*idx2';
    else
        final_idx2=[1:n2]';
        idx1=1:n1;
        final_idx1=S*idx1';
    end
    corr=[final_idx1,final_idx2];
    if n1<=n2
        final_idx1_2 = [1:n1]';
        idx2=1:n2;
        final_idx2_2 = S2'*idx2';
    else
        final_idx2_2=[1:n2]';
        idx1=1:n1;
        final_idx1_2=S2*idx1';
    end
    corr2=[final_idx1_2,final_idx2_2]; 
    errors = zeros(size(corr,1), 1);
    
    for m=1:size(corr,1)
        
        if (strcmp(track, 'low resolution/'))
            gt_match = gt(gt(:,1) == corr(m,1), 2);
            match = corr(m,2);
            
            if ~isempty(gt_match)&& match>0
                % using the second graph
                errors(m) = dist(gt_match, match); % TODO include your geodesics here
            else
                errors(m) = -2;
            end
        else
            errors(m) = -1;
        end
        
    end
    diameters = sqrt(sum(calc_tri_areas(N)));
    errors = errors / diameters;
    for m=1:length(thresholds)
        curve1(k,m) = 100*sum(errors <= thresholds(m)) / length(errors);
    end
    
    errors2 = zeros(size(corr2,1), 1);
    
    for m=1:size(corr2,1)
        
        if (strcmp(track, 'low resolution/'))
            gt_match = gt(gt(:,1) == corr2(m,1), 2);
            match = corr2(m,2);
            
            if ~isempty(gt_match)&& match>0
                % using the second graph
                errors2(m) = dist(gt_match, match); % TODO include your geodesics here
            else
                errors2(m) = -2;
            end
        else
            errors2(m) = -1;
        end
        
    end
    diameters = sqrt(sum(calc_tri_areas(N)));
    errors2 = errors2 / diameters;
    for m=1:length(thresholds)
        curve2(k,m) = 100*sum(errors2 <= thresholds(m)) / length(errors2);
    end
end