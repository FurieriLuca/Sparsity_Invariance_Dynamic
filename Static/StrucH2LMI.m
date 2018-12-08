function [K,J,Jdiag] = StrucH2LMI(A,B1,B2,Gp,Q,R,SP)
% Structured Optimal control over directed graphs: SDP relaxation via block
% diagnal Lyapunov function
% Input data: graph Gp, Gc; Dynamic matrices: A, B1, B2, --> cell format
%             Performance Index: Q1,Q2 --> penalize on state, Q1: absolute, Q2: relative
% Outpute data: Jdiag, performance, K, conresponding controller

epsilon = 0.000001;

%% Obtain dynamics matrices
[N,~] = size(Gp);               % Number of nodes in the graph
[n,m] = size(B2{1});            % dimensions of input and state in each node 

[Amat, Bmat1, Bmat2] = NetStateModel(A,B1,B2,Gp); 

%% solution via Yalmip
% variables
X = sdpvar(n);               %% block diagonal X 
for i = 2:N
    X = blkdiag(X,sdpvar(n));
end

%{
Z = sdpvar(N*m,N*n);        %% Matrix Z has sparsity pattern in Gc
for i = 1:N                     
    for j = 1:N
        if Gc(i,j) == 0 &&  i ~= j
             Z((i-1)*m+1:i*m,(j-1)*n+1:j*n) = zeros(m,n);
        end
    end
end
%}
Z = sdpvar(N*m,N*n);        %% Matrix Z has sparsity pattern in Gc
for i = 1:N*m                     
    for j = 1:N*n
        if SP(i,j)==0
             Z(i,j) = 0;
        end
    end
end

% constraint
Y = sdpvar(N*m);
Const = [X-epsilon*eye(n*N) >=0, ...
    (Amat*X-Bmat2*Z)+(Amat*X-Bmat2*Z)'+Bmat1*Bmat1' + epsilon* eye(n*N) <= 0,...
    [Y Z; Z' X]>=0];

% cost function
Obj = trace(Q*X)+trace(R*Y);

ops = sdpsettings('solver','sedumi');
Info = optimize(Const,Obj,ops);
% solution

X1 = value(X);Z1 = value(Z);Y1 = value(Y);

K = Z1*X1^(-1);  % controller
Jdiag = sqrt(trace(Q*X1)+trace(R*Y1));

%% H2 performance using Lyapunov equation
P = lyap((Amat-Bmat2*K)',Q+K'*R*K);
J = sqrt(trace(P*(Bmat1*Bmat1')));

end

