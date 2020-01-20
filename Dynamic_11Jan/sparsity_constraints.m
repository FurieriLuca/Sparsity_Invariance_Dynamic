
s = sym('s');                                        % Redefines "s" as symbolic variable instead of a tf variable, which would be wrongly interpreted by matlab later
assume(s,'real');
Gs = [1/(s+1) 0 0 0 0;                               % Redefines the plant in terms of symbolic "s"
    1/(s+1) 1/(s-1) 0 0 0;
    1/(s+1) 1/(s-1) 1/(s+1) 0 0;
    1/(s+1) 1/(s-1) 1/(s+1) 1/(s+1) 0;
    1/(s+1) 1/(s-1) 1/(s+1) 1/(s+1) 1/(s-1)];
P21s   = [Gs eye(n)];

I = eye(n);
% Defines CQ and DQ both as symbolic variables and sdpvar variables
CQs = sym('CQ',[m size(P21s,2)*N]);
DQs = sym('DQ',[m size(P21s,2)]);
CQv = sdpvar(m,size(P21s,2)*N);          % decision variables
DQv = sdpvar(m,size(P21s,2));

P21s_inv=simplify(P21s'*inv(P21s*P21s'));

%Force Y in RHinfinity
for(i=1:m)
    for(j=1:size(P21s,2))
        Y(i,j)=CQs(i,[(j-1)*N+1:j*N])*inv(s*eye(N)-AiQ)*BiQ+DQs(i,j);
    end
end

Constraints = [];

%% Constraint Type 1: Y_tilda(s)=Y(s)P21_inv(s) \in Sparse(T) in terms of CQv and DQv
%{
fprintf('Step 1: Encoding the constraint Y_tilda(s) in T ...')
Y_tilda=Y*P21s_inv;
for i = 1:m
    for j = 1:n
        fprintf('   Percentage %6.4f \n', 100*(2*n*(i-1)+j)/m/n/2 );
        
        if Tbin(i,j) == 0     % Whenever we need Y_tilda(i,j) = 0 ....
            [num,~] = numden(Y_tilda(i,j));
            cc      = coeffs(num,s);                                  % All elements of this vector must be 0....
            A_eq    = equationsToMatrix(cc,[vec(CQs);vec(DQs)]);      % Express system of equations in matrix form in terms of the vectorized versions of CQs and DQs
            A_eqs   = double(A_eq);    %A_eqs is the same as A_eq, for computation with sdpvars
            
            Constraints = [Constraints, A_eqs*[vec(CQv);vec(DQv)] == 0]; % Add the constraints in terms of the sdpvars CQv and DQv, by using A_eqs computed with symbolics
        end
    end
end
fprintf('Done \n')
fprintf('Step 2: Encoding the constraint X(s) in R ...')
%}

%% Constraint Type 2: G(s)Y(s) \in Sparse(R) in terms of CQv and DQv
GY_tilda=Gs*Y*P21s_inv;
if QI == 0        %This cycle is useless if QI (redundant constraints). Hence, we skip it in this case.
    %% precomputation
    %%
    for i = 1:n
        for j = 1:n
            fprintf('   Percentage %6.4f \n', 100*(n*(i-1)+j)/n/n );
            if Rbin(i,j) == 0     % Whenever we need GY(i,j) = 0 ....
                [num,~] = numden(GY_tilda(i,j));
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
fprintf('Step 3: Encoding stability of X ...\n')

%% Constraint Type 3: P21+GY is stable
% Defines CQ2 and DQ2 both as symbolic variables and sdpvar variables
CQ2s = sym('CQ2',[m size(P21s,2)*N]);
DQ2s = sym('DQ2',[m size(P21s,2)]);
CQ2v = sdpvar(m,size(P21s,2)*N);          % decision variables
DQ2v = sdpvar(m,size(P21s,2));


%Force Y_aux in RHinfinity
for(i=1:m)
    for(j=1:size(P21s,2))
        Y_aux(i,j)=CQ2s(i,[(j-1)*N+1:j*N])*inv(s*eye(N)-AiQ)*BiQ+DQ2s(i,j);
    end
end
equation=P21s+Gs*Y-Y_aux; %will set this to 0, so P21+GY is stable
for i = 1:m
    for j = 1:2*n
        fprintf('   Percentage %6.4f \n', 100*(2*n*(i-1)+j)/m/n/2 );
        [num,~] = numden(equation(i,j));
        cc      = coeffs(num,s);                                  % All elements of this vector must be 0....
        [A_eq,b_eq]    = equationsToMatrix(cc,[vec(CQs);vec(DQs);vec(CQ2s);vec(DQ2s)]);      % Express system of equations in matrix form in terms of the vectorized versions of CQs and DQs
        A_eqs   = double(A_eq);    %A_eqs is the same as A_eq, for computation with sdpvars
        b_eqs    = double(b_eq);
        Constraints = [Constraints, A_eqs*[vec(CQv);vec(DQv);vec(CQ2v);vec(DQ2v)]== b_eqs]; % Add the constraints in terms of the sdpvars CQv and DQv, by using A_eqs computed with symbolics
    end
end
fprintf('Done \n')

