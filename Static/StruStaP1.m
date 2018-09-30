function [K] = StruStaP1(A,B1,B2,Gp,Gc)
% Solution to the structured stablization P1

%% Obtain dynamics matrices
    [N,~] = size(Gp);               % Number of nodes in the graph
    [n,m] = size(B2{1});            % dimensions of input and state in each node 
    [Amat,~, Bmat2] = NetStateModel(A,B1,B2,Gp); %% Reformulation of data --> the whole state space model

    epsilon = 1.0e-3;
    %% variables
    X = sdpvar(n);               %% block diagonal X 
    for i = 2:N
        X = blkdiag(X,sdpvar(n));
    end
    
    Z = sdpvar(N*m,N*n);        %% Matrix Z has sparsity pattern in Gc
    for i = 1:N
        for j = 1:N
            if Gc(i,j) == 0 && i~=j
                Z((i-1)*m+1:i*m,(j-1)*n+1:j*n) = 0;
            end
        end
    end

    %% constraint
    Const = [X - epsilon*eye(n*N)>= 0, (Amat*X-Bmat2*Z)+(Amat*X-Bmat2*Z)' + epsilon* eye(n*N) <= 0];
    
    %% cost function
    Obj = 0;

    %% solution via sedumi
    ops = sdpsettings('solver','sedumi');
    Info = optimize(Const,Obj,ops);
    X = value(X);Z = value(Z);
    K = Z*X^(-1);
end

