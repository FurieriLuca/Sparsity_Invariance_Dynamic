function [K,J,X1,Y1,Z1,J_restriction] = StrucH2LMI_new_Gamma(A,B1,B2,Q,R,T,Rstruct,Gamma)
% Structured Optimal control over directed graphs: SDP relaxation via block
% diagnal Lyapunov function
% Input data: graph Gp, Gc; Dynamic matrices: A, B1, B2, --> cell format
%             Performance Index: Q1,Q2 --> penalize on state, Q1: absolute, Q2: relative
% Outpute data: Jdiag, performance, K, conresponding controller

epsilon = 1e-6;

n=size(A,1);
m=size(B2,2);

Const=[];
%% solution via Yalmip
% variables
X = sdpvar(n,n);        %% Matrix Z has sparsity pattern in Gc


Y = sdpvar(m,n,'full');        %% Matrix Z has sparsity pattern in Gc

Z = sdpvar(m,m);
% constraint

Const = [norm((X*Gamma).*(1-Rstruct),2)<= 1.0e-8, norm((Y*Gamma).*(1-T),2)<= 1.0e-8, X-epsilon*eye(n) >=0, ...
    (A*X+B2*Y)+(A*X+B2*Y)'+B1*B1' + epsilon* eye(n) <= 0,...
    [Z Y; Y' X]>=0];

% cost function
Obj = trace(Q*X)+trace(R*Z);

ops = sdpsettings('solver','sedumi');
%ops = sdpsettings('solver','sdpt3','verbose',1);

Info = optimize(Const,Obj,ops);
% solution

J_restriction  = value(Obj);%trace(Q*X1)+trace(R*Z1);


X1 = value(X);Z1 = value(Z);Y1 = value(Y);

K = Y1*X1^(-1);  % controller

%% H2 performance using Lyapunov equation
P = lyap((A+B2*K)',Q+K'*R*K);
J = sqrt(trace(P*(B1*B1')));

end

