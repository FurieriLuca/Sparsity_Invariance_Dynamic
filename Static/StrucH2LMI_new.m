function [K,J,Jdiag,r] = StrucH2LMI_new(A,B1,B2,Gp,Q,R,SP)
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

R_struct = generate_SXlessS(SP);                         % This is the MSI matrix, non-symmetric
R_struct = antisymmetrize_bin(R_struct)               % This takes the largest symmetrix component

r = degree_separation(R_struct);

%% solution via Yalmip
% variables
X = sdpvar(N*n,N*n);               %% block diagonal X 

for i = 1:(N*n)                             %%%% X has the structure of R    
    for j = 1:(N*n)
        if R_struct(i,j) == 0
             X(i,j) = 0;
        end
    end
end

%B=ones(n,n)      %% block diagonal X 
%for i = 2:N
%end


Z = sdpvar(N*m,N*n);        %%  Z has the sparsity of T = SP
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

