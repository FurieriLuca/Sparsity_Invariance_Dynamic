% created 2010-May-05
% solve the linear equation for feedback gain F in the alternating
% iterative method for the minimization of augmented Lagrangian
% Here, if R is a diagonal matrix, we exploit the block diagonal structure
% of the corresponding BigA = kron(R,C2*L*C2').
% If R is not a diagonal matrix, the code inverts the matrix BigA, which is
% clearly not efficient. In this case, conjugate gradient method to solve
% this linear equation is recommended.

% note that the original big A matrix for the linear equation is 
% BigA = kron(C2*L*C2',R), which is not a block-diagonal matrix. Therefore,
% we solve the linear equation for transpose F.

function F = SH2_AugLag_ALT_solveF(P,L,V,B2,C2,R,ISc,c)


    [m,n] = size(V);

% if R is diagonal matrix, then the left-hand-side is a block diagonal
% matrix, i.e., kron(R,C2*L*C2') is block diagonal
% then we can solve the equation block by block

if norm(R - R.*eye(size(R))) == 0 
    
    diagR = diag(R); 
    
    % the number of blocks 
    
    r = length(diagR);
    
    % the size of each block is p-by-p
    
    C2LC2 = C2*L*C2';
    p = length(C2LC2);   
        
    vecISct = vec(ISc');
    
    % vector form of the right hand side
    
    vecb = vec(2*C2*L*P*B2 - V');
    
    % declare variable for vector form F transpose
    
    vecFt = zeros(m*n,1);

    % solve all r blocks one by one
    
    for i=1:r
        vecFt((i-1)*p+1:i*p) = (2*diagR(i)*C2LC2 + ...
            c*diag(vecISct((i-1)*p+1:i*p)))\vecb((i-1)*p+1:i*p);
    end

%     This two lines are to check if the solution matches with the result
%     from direction computation
%     vecFt1 = (2*kron(R,C2LC2) + c*diag(vecISct))\vecb;
%     norm(vecFt-vecFt1)

    % put the vector form into matrix form

    Ft = reshape(vecFt,n,m);
    
    % take the transpose to obtain the real 
    
    F = Ft';
else
    % if R is not diagonal matrix, use direct method to solve the
    % large-scale system, but we should try other method to reduce the
    % computational complexity
    
    Abig = 2*kron(C2*L*C2',R) + c*diag(vec(ISc));
    bbig = vec(2*B2'*P*L*C2' - V);
    F    = reshape(Abig\bbig,m,n);
    
end

% check the solution of F

Res     = 2*(R*F*C2-B2'*P)*L*C2' + V + c*(F.*ISc);
normRes = norm(Res,'fro');

% check if normRes is a number, because it is possible there is no solution
% of the linear equation, i.e., BigA is not invertible

if isnan(normRes)==1   
    error('No solution of the linear equation for F!')
end

% check the solution

if normRes > 1.e-5
    disp('check the solution of feedback gain F from solveF!!')
end

