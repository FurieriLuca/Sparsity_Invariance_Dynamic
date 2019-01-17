
s = sym('s');                                        % Redefines "s" as symbolic variable instead of a tf variable, which would be wrongly interpreted by matlab later
Gs = [1/(s+1) 0 0 0 0;                               % Redefines the plant, Knom and the factorization in terms of symbolic "s"
    1/(s+1) 1/(s-1) 0 0 0;
    1/(s+1) 1/(s-1) 1/(s+1) 0 0;
    1/(s+1) 1/(s-1) 1/(s+1) 1/(s+1) 0;
    1/(s+1) 1/(s-1) 1/(s+1) 1/(s+1) 1/(s-1)];
Knoms = [0 0 0 0 0;0 -6/(s+3) 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 0 0 -6/(s+3)];
M_tilda  = inv(eye(n)-Gs*Knoms);
M  = -(eye(m)-Knoms*Gs);
N_tilda  = Gs*inv(eye(m)-Knoms*Gs);
N  = -Gs*inv(eye(m)-Knoms*Gs);
V  = -eye(m);
U  =  Knoms;
V_tilda  = eye(n);
U_tilda  = -Knoms; 

I = eye(n);
% Defines CQ and DQ both as symbolic variables and sdpvar variables
CQs = sym('CQ',[m n*order]);
DQs = sym('DQ',[m n]);
CQv = sdpvar(m,n*order);          % decision variables
DQv = sdpvar(m,n);

Constraints = [];

for(i=1:m)                      % express Y(s) in terms of decision variables (Remark 5 of [2])
    for(j=1:n)
        Y(i,j) =  CQs(i,(j-1)*order+1:j*order)*inv(s*eye(order)-AiQ)*BiQ+DQs(i,j);
    end
end

U_tildaPlusMYTimesM_tilda=(U_tilda+M*Y)*M_tilda;
V_tildaMinusNYTimesM_tilda=(V_tilda-N*Y)*M_tilda;

%% Constraint Type 1: Y(s) \in Sparse(T) in terms of CQv and DQv (same as [2])
fprintf('Step 1: Encoding the constraint Y(s) in T ...')
for i = 1:m
    for j = 1:n
        if Tbin(i,j) == 0     
            %Uij = ;     % Performs the symbolic operation U_tilda*M_tilda+M*w GY(i,j)=\sum_l G(i,l)Y(l,j)
            [num,~] = numden(U_tildaPlusMYTimesM_tilda(i,j));
            cc      = coeffs(num,s);                                  % All elements of this vector must be 0....
            A_eq    = equationsToMatrix(cc,[vec(CQs);vec(DQs)]);      % Express system of equations in matrix form in terms of the vectorized versions of CQs and DQs
            A_eqs   = double(A_eq);    %A_eqs is the same as A_eq, for computation with sdpvars
            
            %for(f=1:size(A_equations,1)) %ALTERNATIVE FORM
            %    for(g=1:n*N)
            %Constraints=[Constraints, A_equations_normal(f,[(g-1)*m+1:g*m])*CQv(:,g)==0];
            %    end
            %end
            
            Constraints = [Constraints, A_eqs*[vec(CQv);vec(DQv)] == 0]; % Add the constraints in terms of the sdpvars CQv and DQv, by using A_eqs computed with symbolics
        end
    end
end
fprintf('Done \n')

%% Constraint Type 2: G(s)Y(s) \in Sparse(R) in terms of CQv and DQv
fprintf('Step 2: Encoding the constraint GY(s) in R ...\n')
if QI == 0        %This cycle is useless if QI (redundant constraints). Hence, we skip it in this case.
    %% precomputation
    Gi = (s*eye(order)-AiQ)\BiQ;
    %%
    for i = 1:n
        for j = 1:n
            fprintf('   Percentage %6.4f \n', (n*(i-1)+j)/n/n );
            if Rbin(i,j) == 0     % Whenever we need GY(i,j) = 0 ....
                 %Uij = Gs(i,:)* (CQs(:,(j-1)*N+1:j*N) * Gi + DQs(:,j));     % Performs the symbolic matrix product GY(i,j)=\sum_l G(i,l)Y(l,j)
                [num,~] = numden(V_tildaMinusNYTimesM_tilda(i,j));
                cc      = coeffs(num,s);                                  % All elements of this vector must be 0....
                A_eq    = equationsToMatrix(cc,[vec(CQs);vec(DQs)]);      % Express system of equations in matrix form in terms of the vectorized versions of CQs and DQs
                A_eqs   = double(A_eq);    %A_eqs is the same as A_eq, for computation with sdpvars
                
                %for(f=1:size(A_equations,1)) %ALTERNATIVE FORM
                %    for(g=1:n*N)
                %Constraints=[Constraints, A_equations_normal(f,[(g-1)*m+1:g*m])*CQv(:,g)==0];
                %    end
                %end
                
                Constraints = [Constraints, A_eqs*[vec(CQv);vec(DQv)] == 0]; % Add the constraints in terms of the sdpvars CQv and DQv, by using A_eqs computed with symbolics
            end
        end
    end
end
fprintf('Encoding the constraint GY(s) in R ...Done\n')