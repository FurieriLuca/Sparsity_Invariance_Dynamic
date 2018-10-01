function [Kopt,Jopt,Iter] = StrucH2_Gradient(A,B1,B2,Gp,Gc,Q,R,K0)
    % the code solves the structured H2 problem using gradient method

    
    %% Converting the model into state space
    [N,~] = size(Gp);               % Number of nodes in the graph
    [n,m] = size(B2{1});            % dimensions of input and state in each node 

    [Amat, Bmat1, Bmat2] = NetStateModel(A,B1,B2,Gp); 
    
    SP = kron(Gc+eye(N),ones(m,n));  %% sparsity patten
    
    %% Initial Stablization Structured controller    
    %K0 = StruStaP1(A,B1,B2,Gp,Gc);
    %[K0,~,~] = StrucH2LMI(A,B1,B2,Gp,Gc,Q,R);
    J0 = trace(lyap((Amat - Bmat2*K0)',Q + K0'*R*K0)*(Bmat1*Bmat1'));
    P = lyap((Amat - Bmat2*K0)',Q + K0'*R*K0); 
    L = lyap((Amat - Bmat2*K0),Bmat1*Bmat1'); 
    gradK0  = 2*(R*K0*L - Bmat2'*P*L).*SP;    %% initial gradient
    NgradK0 = norm(gradK0,'fro'); 
    
    %% Parameters in iteration  
    alpha = 0.01;     % backtrapping line search 
    beta  = 0.5;
    tol   = 0.00001;    % tolerance of norm of gradient direction
    MaxIter = 500;    % maximum number of gradient steps
    Disp = 10;
    
    K = K0; J = J0;     %% starting point 
    disp('Iter       ngradK           J') 
    for Iter = 1:MaxIter  %% iterations of gradient projection

        % compute the gradient projection
        P = lyap((Amat - Bmat2*K)',Q + K'*R*K); 
        L = lyap((Amat - Bmat2*K),Bmat1*Bmat1'); 
        gradK = 2*(R*K*L - Bmat2'*P*L);    %% gradient
        PgradK = gradK.*SP;                %% projected gradient
        ngradK = norm(PgradK,'fro');
        
        if mod(Iter,Disp) == 0 || Iter == 1
            disp([num2str(Iter),'         ',num2str(ngradK,'%6.4E'),'     ',num2str(J,'%6.4E')])
        end
        
        % update according to Armijo rule
        SteSize = 1;
        Ktemp = K - SteSize*PgradK;  %% update a step
        maxEigAc = max(real(eig(Amat - Bmat2*Ktemp)));
        Jtemp    = trace(lyap((Amat - Bmat2*Ktemp)',Q + Ktemp'*R*Ktemp)*(Bmat1*Bmat1'));
        
        while maxEigAc >= 0 || J - Jtemp < alpha*trace(gradK'*(SteSize*PgradK))
            SteSize = beta*SteSize;
            if SteSize < 1.e-19
                disp('Gradient method gets stuck with very small stepsize')
            end
            Ktemp = K - SteSize*PgradK;
            maxEigAc = max(real(eig(Amat-Bmat2*Ktemp)));
            Jtemp = trace(lyap((Amat-Bmat2*Ktemp)', Q + Ktemp'*R*Ktemp)*(Bmat1*Bmat1'));
            break
        end
        
        % update the current step K
        K = Ktemp;
        J = Jtemp;
        
        % stop the algorithm if the norm of gradient is small enough        
        if ngradK < tol %ngradK / NgradK0 < tol
            break;
        end
    end
    
    Kopt = K;
    Jopt = sqrt(J);
        
end

