function [A, B, P_rnd] = generate_wig_2(n,sigma,returnid)
    %A=normrnd(0,var/n,n);
    %A=mvnrnd(zeros(n,1),eye(n)/n,n);
    A=randn(n)*(1/sqrt(n));
    %A=tril(A,-1)+tril(A,-1)'+diag(normrnd(0,2*var/n,n,1));
    %A=tril(A,-1)+tril(A,-1)'+diag(mvnrnd(zeros(n,1),(2/n)*eye(n),1));
    aux_diag=randn(n,1)*(sqrt(2)/sqrt(n));
    A=tril(A,-1)+tril(A,-1)'+diag(aux_diag);
    
    %B=normrnd(0,var/n,n);
    %B=mvnrnd(zeros(n,1),eye(n)/n,n);
    B=randn(n)*(1/sqrt(n));
    %B=tril(B,-1)+tril(B,-1)'+diag(normrnd(0,2*var/n,n,1));
    %B=tril(B,-1)+tril(B,-1)'+diag(mvnrnd(zeros(n,1),(2/n)*eye(n),1));
    aux_diag=randn(n,1)*(sqrt(2)/sqrt(n));
    B=tril(B,-1)+tril(B,-1)'+diag(aux_diag);
    B=A+sigma*B;
    P_rnd = eye(n);
    if returnid== 0
        P_rnd = P_rnd(:, randperm(n));
        B = P_rnd * B * P_rnd';
    end
    