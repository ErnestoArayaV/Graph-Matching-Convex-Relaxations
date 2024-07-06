function [P,diagdom,prop] = matching_PGD(A,B,max_iter,str)
%projected gradient descent alg. for graph matching
%
if strcmp(str,'dynamic_steps')
    [P,diagdom,prop] = matching_pgd_dynsteps_loop(A,B,max_iter);
    return 
elseif strcmp(str,'fixed_steps')
    [P,diagdom,prop] = matching_pgd_fixsteps_loop(A,B,max_iter);
    return
elseif strcmp(str,'heuristic')
    [P,diagdom,prop] = matching_pgd_heuristic_loop(A,B,max_iter);
end
end