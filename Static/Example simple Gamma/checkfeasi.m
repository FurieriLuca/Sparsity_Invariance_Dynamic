function [flag,Cost] = checkfeasi(X,Y,Z,A,B1, B2, Q,R, Gamma,T,Rs)

epsilon = 1e-3;
n = size(A);

% sparsity pattern
    flag = 0;
    if norm((X*Gamma).*(1-Rs),2)> 1.0e-8
        flag = flag + 1;
    end

    if norm((Y*Gamma).*(1-T),2)> 1.0e-8
        flag = flag + 1;
    end
    
    temp = [Z Y;Y' X];
    
    if min(eig(temp)) < -1.0e-6
        flag = flag + 1;
    end
    temp =  (A*X+B2*Y)+(A*X+B2*Y)'+B1*B1' + epsilon* eye(n);
    if max(eig(temp))  >= 1.0e-6
        flag = flag + 1;
    end
    

    if flag == 0
        Cost = trace(Q*X) + trace(R*Z);
    else
        Cost = inf;
    end


end

