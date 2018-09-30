function [ Amat, Bmat1, Bmat2] = NetStateModel(A,B1,B2,Gp)
% Convert the data into state space form
%   Detailed explanation goes here

    [N,~] = size(Gp);     % Number of nodes in the graph
    [n,m] = size(B1{1});      % dimensions of input and state in each node 

    Amat = zeros(N*n); Bmat1 = zeros(N*n,N*m);
    for i = 1:N
        for j = 1:N
            if i == j
                Amat((i-1)*n+1:i*n,(j-1)*n+1:j*n) = A{i,j};
                Bmat1((i-1)*n+1:i*n,(j-1)*m+1:j*m) = B1{i};
                Bmat2((i-1)*n+1:i*n,(j-1)*m+1:j*m) = B2{i};
            elseif Gp(i,j) == 1
                Amat((i-1)*n+1:i*n,(j-1)*n+1:j*n) = A{i,j};
            end
        end
    end
end

