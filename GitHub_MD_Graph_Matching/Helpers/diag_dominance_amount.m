function quant=diag_dominance_amount(X)
n= size(X,1);
[~,ind_X]=max(X,[],2);
quant=n-size(unique(ind_X),1);