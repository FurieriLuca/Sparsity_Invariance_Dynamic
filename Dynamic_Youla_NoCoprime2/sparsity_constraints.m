
s = sym('s');                                        % Redefines "s" as symbolic variable instead of a tf variable, which would be wrongly interpreted by matlab later
Gs = [1/(s+1) 0 0 0 0;                               % Redefines the plant, Knom and the factorization in terms of symbolic "s"
        1/(s+1) 1/(s-1) 0 0 0;
        1/(s+1) 1/(s-1) 1/(s+1) 0 0;
        1/(s+1) 1/(s-1) 1/(s+1) 1/(s+1) 0;
        1/(s+1) 1/(s-1) 1/(s+1) 1/(s+1) 1/(s-1)];


I = eye(n);
% Defines CQ and DQ both as symbolic variables and sdpvar variables
CQs = sym('CQ',[m n*order]);
DQs = sym('DQ',[m n]);
CQv = sdpvar(m,n*order);          % decision variables
DQv = sdpvar(m,n);

CQs2 = sym('CQ2',[n n*order]);
DQs2 = sym('DQ2',[n n]);
CQv2 = sdpvar(n,n*order);          % decision variables
DQv2 = sdpvar(n,n);

Constraints = [];

for(i=1:m)                      % express Y(s) in terms of decision variables (Remark 5 of [2])
        for(j=1:n)
               Y(i,j) =  CQs(i,(j-1)*order+1:j*order)*inv(s*eye(order)-AiQ)*BiQ+DQs(i,j);
        end
end

for(i=1:n)                      % express X(s) in terms of decision variables (Remark 5 of [2])
        for(j=1:n)
                X(i,j) =  CQs2(i,(j-1)*order+1:j*order)*inv(s*eye(order)-AiQx)*BiQ+DQs2(i,j);
        end
end


%% Constraint Type 1: Y(s) \in Sparse(T) and X(s) \in Sparse(R) in terms of CQv, DQv and CQv2, DQv2 (similar to [2])
fprintf('Step 1: Encoding the constraint Y(s) in T ...')
for i = 1:m
        for j = 1:n
                if Tbin(i,j) == 0 % main cycle
                        Constraints = [Constraints, CQv(i,[(j-1)*order+1:j*order]) == 0, DQv(i,j) == 0];
                end
        end
end
fprintf('Done \n')

fprintf('Step 2: Encoding the constraint X(s) in R ...\n')
for i = 1:n
        for j = 1:n
                if Rbin(i,j) == 0 % main cycle
                        Constraints = [Constraints, CQv2(i,[(j-1)*order+1:j*order]) == 0, DQv2(i,j) == 0];
                end
        end
end
fprintf('Done \n')

%% Constraint Type21: Achievability of the corresponding controller
fprintf('Step 3: Encoding the constraint X(s)=I-G(s)*Y(s) ...\n')
achievability=X-I+Gs*Y;
for i = 1:n
        for j = 1:n
                fprintf('   Percentage %6.4f \n', 100*(n*(i-1)+j)/n/n );
                [num,~] = numden(achievability(i,j));
                cc      = coeffs(num,s);                                  % All elements of this vector must be 0....
                [A_eq,b_eq]    = equationsToMatrix(cc,[vec(CQs);vec(DQs);vec(CQs2);vec(DQs2)]);      % Express system of equations in matrix form in terms of the vectorized versions of CQs and DQs
                A_eqs   = double(A_eq);    %A_eqs is the same as A_eq, for computation with sdpvars
                b_eqs    = double(b_eq);
                Constraints = [Constraints, A_eqs*[vec(CQv);vec(DQv);vec(CQv2);vec(DQv2)]== b_eqs]; % Add the constraints in terms of the sdpvars CQv and DQv, by using A_eqs computed with symbolics
        end
end
fprintf('Done\n')
