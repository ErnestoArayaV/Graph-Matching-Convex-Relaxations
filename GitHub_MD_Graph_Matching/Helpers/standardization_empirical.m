%standardization of graphs ott be matched. Used in autonomous systems data.
%See our paper for details.
function mtx_stan = standardization_empirical(mtx)
mtx_stan = mtx- mean(mtx,'all');
mtx_stan = mtx_stan/std(mtx,0,'all');
