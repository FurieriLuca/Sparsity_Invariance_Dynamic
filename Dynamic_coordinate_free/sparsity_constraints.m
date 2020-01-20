

Gnom=G*inv(eye(n)-Knom*G);

I = eye(n);
% Defines CQ and DQ both as symbolic variables and sdpvar variables
CQs11 = sym('CQ11',[n n*N]);
DQs11 = sym('DQ11',[n n]);
CQv11 = sdpvar(n,n*N);
DQv11 = sdpvar(n,n);

CQs12 = sym('CQ12',[n m*N]);
DQs12 = sym('DQ12',[n m]);
CQv12 = sdpvar(n,m*N);
DQv12 = sdpvar(n,m);

CQs21 = sym('CQ21',[m n*N]);
DQs21 = sym('DQ21',[m n]);
CQv21 = sdpvar(m,n*N);
DQv21 = sdpvar(m,n);

CQs22 = sym('CQ22',[m m*N]);
DQs22 = sym('DQ22',[m m]);
CQv22 = sdpvar(m,m*N);
DQv22 = sdpvar(m,m);

Constraints = [];

%%CREATES Q11,Q12,Q21,Q22


for i = 1:n
    for j = 1:n
        Q11(i,j) = CQs11(i,(j-1)*N+1:j*N)*inv(s*eye(size(AiQ,1))-AiQ)*BiQ + DQs11(i,j);
    end
end

for i = 1:n
    for j = 1:m
        Q12(i,j) = CQs12(i,(j-1)*N+1:j*N)*inv(s*eye(size(AiQ,1))-AiQ)*BiQ + DQs12(i,j);
    end
end

for i = 1:m
    for j = 1:n
        Q21(i,j) = CQs21(i,(j-1)*N+1:j*N)*inv(s*eye(size(AiQ,1))-AiQ)*BiQ + DQs21(i,j);
    end
end

for i = 1:m
    for j = 1:m
        Q22(i,j) = CQs22(i,(j-1)*N+1:j*N)*inv(s*eye(size(AiQ,1))-AiQ)*BiQ + DQs22(i,j);
    end
end

expression = simplifyFraction(Knom*Q11+Knom*Q12*Knom+Q21+Q22*Knom);
[num,den]=numden(expression);
%% Constraint Type 1: Y(s) \in Sparse(T) in terms of CQv and DQv (same as [2])
fprintf('Step 1: Encoding the constraint Q(s) in S ...')
for i = 1:m
    for j = 1:n
        if Tbin(i,j) == 0 
            cc      = coeffs(num(i,j),s);                                  
            [A_eq,b_eq]    = equationsToMatrix(cc,[vec(CQs11);vec(DQs11);vec(CQs12);vec(DQs12);vec(CQs21);vec(DQs21);vec(CQs22);vec(DQs22)]);      
            A_eqs   = double(A_eq);    
            b_eqs = double(b_eq);
            Constraints = [Constraints, A_eqs*[vec(CQv11);vec(DQv11);vec(CQv12);vec(DQv12);vec(CQv21);vec(DQv21);vec(CQv22);vec(DQv22)] == b_eqs]; % Add the constraints in terms of the sdpvars CQv and DQv, by usin
        end
    end
end
fprintf('Done \n')

