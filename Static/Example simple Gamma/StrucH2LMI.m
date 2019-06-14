function [K,J] = StrucH2LMI(A,B1,B2,Q,R,S)
% Structured Optimal control over directed graphs: SDP relaxation via block
% diagnal Lyapunov function
% Input data: graph Gp, Gc; Dynamic matrices: A, B1, B2, --> cell format
%             Performance Index: Q1,Q2 --> penalize on state, Q1: absolute, Q2: relative
% Outpute data: Jdiag, performance, K, conresponding controller

epsilon = 0.00001;

n=size(A,1);
m=size(B2,2);
%% solution via Yalmip
% variables
X = sdpvar(n);               %% block diagonal X 
X=X.*eye(n);

Z = sdpvar(m,n);        %% Matrix Z has sparsity pattern in Gc
for i = 1:m                     
    for j = 1:n
        if S(i,j)==0
             Z(i,j) = 0;
        end
    end
end

% constraint
Y = sdpvar(m);
Const = [X-epsilon*eye(n) >=0, ...
    (A*X-B2*Z)+(A*X-B2*Z)'+B1*B1' + epsilon* eye(n) <= 0,...
    [Y Z; Z' X]>=0];

% cost function
Obj = trace(Q*X)+trace(R*Y);

ops = sdpsettings('solver','sedumi');
Info = optimize(Const,Obj,ops);
% solution

X1 = value(X);Z1 = value(Z);Y1 = value(Y);

K = Z1*X1^(-1);  % controller

%% H2 performance using Lyapunov equation
P = lyap((A-B2*K)',Q+K'*R*K);
J = sqrt(trace(P*(B1*B1')));

end

