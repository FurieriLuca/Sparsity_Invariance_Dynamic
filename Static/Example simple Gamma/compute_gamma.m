function [ Gamma,X,Y] = compute_gamma( A,B1,B2,K,T,Rstruct,Q,R )

n=size(A,1);
m=size(B2,2);
X=sdpvar(n,n);

ops = sdpsettings('solver','sedumi');

epsilon=0.1;
Z=sdpvar(m,m);
optimize([X-epsilon*eye(n)>=0,X*(A+B2*K)'+(A+B2*K)*X+B1*B1'+epsilon*eye(n)<=0,[Z K*X;(K*X)' X]>=0],trace(Q*X)+trace(R*Z),ops); %here I compute X to minimize the restriction cost
X=value(X);
Y=K*X;

Gamma=sdpvar(n,n);


for(i=1:m)
    for(j=1:n)
        if(T(i,j)==1) .  %computes complementary sparsities to use dot product later
            Tc(i,j)=0;
        else
            Tc(i,j)=1;
        end
    end
end
for(i=1:n)
    for(j=1:n)
        if(Rstruct(i,j)==1)
            Rc(i,j)=0;
        else
            Rc(i,j)=1;
        end
    end
end


optimize([(Y*Gamma).*Tc==zeros(m,n),(X*Gamma).*Rc==zeros(n,n),Gamma>=0],1,ops); %Computes Gamma for sparsity, to include K in the restriction
Gamma=value(Gamma);


end

