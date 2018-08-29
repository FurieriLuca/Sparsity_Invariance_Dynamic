function [ QI] = test_QI( Sbin,Delta)

SDS=bin(Sbin*Delta*Sbin);
QI=1;
for(i=1:size(Sbin,1))
    for(j=1:size(Sbin,2))
        if(SDS(i,j)==1 && Sbin(i,j)==0)
            QI=0;
        end
    end
end

end

