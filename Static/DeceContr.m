function Kd = DeceContr(A,B1,B2,Gp,flag)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


N = length(B1);

OutDegree = sum(Gp);

Kd = [];

for i = 1:N
    Ind = find(Gp(i,:) == 1);
    tmp = 0;
    for j = 1:length(Ind)
        tmp = A{i,Ind(j)}*A{i,Ind(j)}'+tmp;
    end
    tmpB = sqrt(tmp);
    
    %% call Yalmip
    
    [n,m] = size(B2{i});
    
    X = sdpvar(n);
    Z = sdpvar(m,n);
    
    epsilon = 1e-4;
%     

if flag == 1
     Const = [[A{i,i}*X - B2{i}*Z + (A{i,i}*X - B2{i}*Z)', tmpB, X; ...
               tmpB', -eye(n), zeros(n);...
               X, zeros(n), -1./sqrt(OutDegree(i))*eye(n)] <= 0];
else
     Const = [[A{i,i}*X - B2{i}*Z + (A{i,i}*X - B2{i}*Z)', tmpB, X; ...
               tmpB', -eye(n), zeros(n);...
               X, zeros(n), -1./(OutDegree(i))*eye(n)] <= 0];
end
    Const = [Const, X - epsilon*eye(n)>= 0];
    Obj = 0;
    
    ops = sdpsettings('solver','sedumi');
    Info = optimize(Const,Obj,ops);
    X = value(X);Z = value(Z);
    K = Z*X^(-1);
    
    Kd = blkdiag(Kd,K);
end
end

