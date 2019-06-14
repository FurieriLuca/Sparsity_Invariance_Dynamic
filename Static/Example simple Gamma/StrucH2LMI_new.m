function [K,J,X1,Y1] = StrucH2LMI_new(A,B1,B2,Q,R,T,Rstruct)
% Structured Optimal control over directed graphs: SDP relaxation via block
% diagnal Lyapunov function
% Input data: graph Gp, Gc; Dynamic matrices: A, B1, B2, --> cell format
%             Performance Index: Q1,Q2 --> penalize on state, Q1: absolute, Q2: relative
% Outpute data: Jdiag, performance, K, conresponding controller

epsilon = 0.0001;

n=size(A,1);
m=size(B2,2);
%% solution via Yalmip
% variables
X = sdpvar(n,n);        %% Matrix X has sparsity pattern in R
for i = 1:n                     
    for j = 1:n
        if Rstruct(i,j)==0
             X(i,j) = 0;
        end
    end
end

Y = sdpvar(m,n);        %% Matrix Y has sparsity pattern in T
for i = 1:m                     
    for j = 1:n
        if T(i,j)==0
             Y(i,j) = 0;
        end
    end
end

% constraint
Z = sdpvar(m);
Const = [X-epsilon*eye(n) >=0, ...
    (A*X+B2*Y)+(A*X+B2*Y)'+B1*B1' + epsilon* eye(n) <= 0,...
    [Z Y; Y' X]>=0];

% cost function
Obj = trace(Q*X)+trace(R*Z);

ops = sdpsettings('solver','sedumi');
Info = optimize(Const,Obj,ops);
% solution

X1 = value(X);Z1 = value(Z);Y1 = value(Y);

K = Y1*X1^(-1);  % controller

%% H2 performance using Lyapunov equation
P = lyap((A+B2*K)',Q+K'*R*K);
J = sqrt(trace(P*(B1*B1')));

end

