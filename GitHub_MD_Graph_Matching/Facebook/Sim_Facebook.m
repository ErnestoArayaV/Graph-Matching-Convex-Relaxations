% Simulations for the MDGM and PGDGM algorithms with the Facebook Data in the Standford project.
%Uses 'Standfor3.mat' Data file. 
%
%%
clear;
addpath(genpath('Helpers'));
addpath(genpath('Algos Graph Matching'));
load Stanford3.mat;
s=0.5:0.1:1;alpha=1; %subsample parameters 
n_subsamples = 15;%number of subsamples 

correct=zeros(3,length(s),n_subsamples);
max_iter = 100;


%% subsampling A
for l=1:n_subsamples
N = 1000;
ind = randperm(size(A,1));
ind = ind(1:N); %random GT permutation
A_sub = A(ind, ind);

for k = 1:length(s)
    
%% subsampling with s and alpha
sample=zeros(N,N);
for i=1:N
    for j=i+1:N
        g=rand<s(k);
        sample(i,j)=g;
        sample(j,i)=g;
    end
end
%W1=A.*sample;
W1 = A_sub.*sample;
t1=rand(1,N)<alpha;
W1(t1==0,:)=0;W1(:,t1==0)=0;
W1=sparse(W1);

t2=rand(1,N)<alpha;
sample=zeros(N,N);
for i=1:N
    for j=i+1:N
        g=rand<s(k);
        sample(i,j)=g;
        sample(j,i)=g;
    end
end
%W2=A.*sample;
W2 = A_sub.*sample;
W2(t2==0,:)=0;W2(:,t2==0)=0;
W2=sparse(W2);


t=t1.*t2;
N1=sum(t);
truth=(1:N).*t;

%% apply algos 

pi_n = matching_robust_spectral(W1,W2,1);
correct(1,k,l)=sum((pi_n==truth).*(pi_n>0))/N1;

pi_n = matching_MD(W1,W2,max_iter,'dynamic_steps');
correct(2,k,l)=sum((pi_n==truth).*(pi_n>0))/N1;

pi_n = matching_PGD(W1,W2,max_iter,'heuristic');
correct(3,k,l)=sum((pi_n==truth).*(pi_n>0))/N1;


end 
end

correct_m = mean(correct,3);
%plots...
figure;hold on;
plot(s,correct_m(2,:),'+--','Color',[0.75 0 0]);%plot md
plot(s,correct_m(1,:),'o--','Color','blue');%plot grampa
plot(s,correct_m(3,:),'^--','Color',[.75 .75 1]);%plot pgd
lgd=legend('EMDGM','Grampa','PGDGM');
lgd.FontSize = 11;
hXLabel = xlabel('$s$','interpreter','latex');
hYLabel = ylabel('Recovery fraction','interpreter','latex');
set([hXLabel, hYLabel], 'FontSize', 20);
