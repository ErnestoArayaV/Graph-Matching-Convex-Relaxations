%computing norms on sample points
path = '/Users/biernack/Documents/MATLAB/MD graph matching/Auto sys/';
load(strcat(path,'/mat_files/auto_sys_data.mat'));

n = 1000;
sampl=10000;
deg = sum(auto_sys_mat{1});
[~, I] = maxk(deg, n);
A = auto_sys_mat{1}(I, I);

norms_inf = zeros(9,1);
norms_sq = zeros(9,1);
for i=1:9
    disp("i" +num2str(i))
    B = auto_sys_mat{i}(I, I);
    largest_norm_inf = 0;
    largest_norm_sq = 0;
    %sample random points 
    for s=1:sampl
        disp(s)
        X = urnd_simp_sparse(n,n,n);
         mataux = A*X-X*B;
        grad = A*mataux-mataux*B;
        grad_norm_inf = max(max(abs(grad)));
        grad_norm_sq = norm(grad,'fro');
        if grad_norm_inf>largest_norm_inf
            largest_norm_inf = grad_norm_inf;
        end
        if grad_norm_sq>largest_norm_sq
            largest_norm_sq = grad_norm_sq;
        end
    end
    norms_inf(i) = largest_norm_inf;
    norms_sq(i) = largest_norm_sq;
end

