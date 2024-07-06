function [diagdom,property]=sufficiency_metrics(X,str)
diagdom = diag_dominance_amount(X); 
[~,property]=sufficient_prop_check(X,str);
