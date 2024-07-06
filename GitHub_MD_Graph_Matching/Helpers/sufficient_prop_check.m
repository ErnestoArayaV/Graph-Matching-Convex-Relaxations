function [mat,nonzeros]=sufficient_prop_check(M,type_string)
n=size(M,1);
mat=zeros(n,n);

if strcmp(type_string,"sum")
    for i=1:n
        for j=1:n
             mat(i,j)=M(i,i)+M(j,j)-M(i,j)-M(j,i);
        end
    end
end

if strcmp(type_string,"max")
    for i=1:n
        for j=1:n
             mat(i,j)=max(M(i,i),M(j,j))-max(M(i,j),M(j,i));
        end
    end
end

if strcmp(type_string,"sum-max")
    for i=1:n
        for j=1:n
             mat(i,j)= 0.5*(M(i,i)+M(j,j))-max(M(i,j),M(j,i));
        end
    end
end

if strcmp(type_string,"prod")
    for i=1:n
        for j=1:n
             mat(i,j)=sqrt(M(i,i)*M(j,j))/M(i,j)-1;
        end
    end
end

aux=mat>0;
nonzeros=(nnz(~aux)-n);
        