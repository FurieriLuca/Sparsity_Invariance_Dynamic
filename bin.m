function [A_bin] = bin( A )

n=size(A,1);
m=size(A,2);
A_bin=zeros(n,m);

for(i=1:n)
    for(j=1:m)
        if(A(i,j)>0.00001 || A(i,j)<-0.00001)
            A_bin(i,j)=1;
        end
    end
end


end

