function [ T] = antisymmetrize_bin( T )


n=size(T,1);
for(i=1:n)
    for(j=1:n)
        if(T(i,j)==0)
            T(j,i)=0;
        end
    end
end

end

