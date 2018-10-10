function [ r] = degree_separation( R )
%Compute how many addends in separation of the Lyapunov function
n = size(R,1);
D=zeros(n,n);

for(i=1:n)
    D(i,i)=sum(R(i,:));
end

L = D - R; %Graph Laplacian
eigs_L=eig(L);

r=0;
for(i=1:n)
        if(eigs_L(i)<0.000001)
                r = r +1;
        end
end

end

