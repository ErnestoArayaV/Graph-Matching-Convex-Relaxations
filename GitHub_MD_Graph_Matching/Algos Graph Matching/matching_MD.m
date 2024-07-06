%Mirror descent alg. for graph matching
%
function [P,diagdom,prop] = matching_MD(A,B,max_iter,str)
if strcmp(str,'dynamic_steps')
    [P,diagdom,prop] = matching_MD_dynsteps_loop(A,B,max_iter);
    return 
elseif strcmp(str,'fixed_steps')
    [P,diagdom,prop] = matching_MD_fixsteps_loop(A,B,max_iter);
end
end

    
% n = size(A,1);
% A2 = A^2;
% B2 = B^2;
% %l_rate = 0.1;
% X = ones(n)/n^2;
% X_best = X;
% largest_grad_norm=0;
% val_truth = 6.6100e-08 ;
% %l_rate = 20;
% %MD loop
% for k=1:max_iter  
%     %print iterations for the impatients 
%     %if mod(k,20) == 0
%     k
%     %l_rate=l_rate/sqrt(k);
%     %end
%     grad = 2*(A2*X+X*B2)-4*A*X*B;
%     grad_norm = max(max(abs(grad)));
%     largest_grad_norm = max([grad_norm,largest_grad_norm]);
%     disp(['grad_norm, ',num2str(grad_norm)])
%     disp(['largest_grad_norm, ',num2str(largest_grad_norm)])
%     
% %    l_rate = 16;
%  % l_rate = (6.6107e-08)
%     if strcmp(str,'dynamic_steps')
%         if grad_norm<0.0000001
%             l_rate = 2*sqrt(2)/(largest_grad_norm*sqrt(k+1));
%         else 
%             l_rate = 2*sqrt(2)/(grad_norm*sqrt(k+1));
% 
%         end 
%         disp(['gamma, ', num2str(l_rate)])
%     end
%     if strcmp(str,'polyak')
%         if grad_norm==0
%             l_rate = 1;
%         else
%             l_rate = (convex_obj(A,B,X)-val_truth)/grad_norm^2;
%         end
%         disp(['gamma',num2str(l_rate)])
%     end
%     %fixed steps
%     if strcmp(str,'fixed_steps')
%         l_rate = sqrt(2*log(n))*largest_grad_norm/sqrt(max_iter+1);
%         disp(['gamma',num2str(l_rate)])
%     end
%     %MD explicit multiplicative updates
%     %grad = B2*X+X*A2-2*B*X*A;
%     egrad = exp(-l_rate*grad);
%     X = X.*egrad;
%     %X = X/sum(sum(X));
%     obj_X = convex_obj(A,B,X);
%     if obj_X<convex_obj(A,B,X_best)
%         X_best = X./sum(sum(X));
%     end
%     disp(['obj, ', num2str(obj_X)])
%     X = X./sum(sum(X));
% end
% %Rounding step
% P = GMWM_alg(X_best',-2000);
% 
