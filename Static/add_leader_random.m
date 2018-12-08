function [leaders,r ] = add_leader_random( leaders,n)

L=size(leaders,2);
assigned=0;

while(assigned==0)
    assigned=1;
    r=randi(n);
    for(i=1:L)
        if(r==leaders(i))
            assigned=0;
        end
    end
end
%%Add leaders randomly
leaders=[leaders,r];



end

