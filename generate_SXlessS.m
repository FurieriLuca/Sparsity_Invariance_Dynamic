function [ X ] = generate_SXlessS( S )
m=size(S,1);
n=size(S,2);
X=ones(n,n);

for(i=1:m)
    for(k=1:n)
        if(S(i,k)==0)
            for(j=1:n)
                if(S(i,j)==1)
                    X(j,k)=0;
                end
            end
        end
    end
end


end

