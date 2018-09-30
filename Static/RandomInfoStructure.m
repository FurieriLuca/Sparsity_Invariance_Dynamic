function [ Gc] = RandomInfoStructure(n)
Gc = zeros(n^2,n^2);

for(i = 1:n^2)
        for(j = 1:n^2)
                if(rand>0.06)
                        Gc(i,j)=1;
                end
        end
end
Gc = Gc + eye(n^2);
end

