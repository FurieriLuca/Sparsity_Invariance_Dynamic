function [] = print_closest_QIs( Sbin,Delta)
%You can use this function to check that with Sbin7 to check that the only
%matrix which is QI and has cardinality greater or equal than 5 is Sbin1.

m=size(Delta,1);
n=size(Delta,2);
NSbin=sum(sum(Sbin));

   count=0;
   count_inside=0;
for(i=2^NSbin:-1:1)
    b=de2bi(i);
    b=[b zeros(1,NSbin-size(b,1))];
    for(j=1:m)
        for(k=1:n)
            if(Sbin(j,k)==1)
                count_inside=count_inside+1;
                T(j,k)=b(count_inside);
            end
        end
    end
    count_inside=0;
    QI=test_QI(T,Delta);
    if(QI==1 && sum(sum(T))>=3)
        count=count+1;
        T
        %if(count==1)
            
       % end
    end
end

end

