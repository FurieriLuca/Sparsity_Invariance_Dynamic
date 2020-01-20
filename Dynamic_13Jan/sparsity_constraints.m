
s = sym('s');                                        % Redefines "s" as symbolic variable instead of a tf variable, which would be wrongly interpreted by matlab later
assume(s,'real');
Gs = [1/(s+1) 0 0 0 0;                               % Redefines the plant in terms of symbolic "s"
    1/(s+1) 1/(s-1) 0 0 0;
    1/(s+1) 1/(s-1) 1/(s+1) 0 0;
    1/(s+1) 1/(s-1) 1/(s+1) 1/(s+1) 0;
    1/(s+1) 1/(s-1) 1/(s+1) 1/(s+1) 1/(s-1)];




CQYs = sym('CQY',[m n*N]);
DQYs = sym('DQY',[m n]);
CQYv = sdpvar(m,n*N);          % decision variables
DQYv = sdpvar(m,n);


CQWs = sym('CQW',[n m*Nw]);
DQWs = sym('DQW',[n m]);
CQWv = sdpvar(n,m*Nw);          % decision variables
DQWv = sdpvar(n,m);



%Force Y in RHinfinity

for(i=1:m)
    for(j=1:n)
        Y(i,j)=CQYs(i,[(j-1)*N+1:j*N])*inv(s*eye(N)-AiQ)*BiQ+DQYs(i,j);
    end
end

for(i=1:n)
    for(j=1:m)
        W(i,j)=CQWs(i,[(j-1)*Nw+1:j*Nw])*inv(s*eye(Nw)-AiQw)*BiQw+DQWs(i,j);
    end
end


Constraints = [];

%% Constraint "Achievability": X=I+GY,W=XG,Z=I+YG

fprintf('Step 1-a: Encoding the constraint X=I+GY is stable \n')
%constraints_X;


fprintf('Step 1-b: Encoding the constraint W=(I+GY)G is stable \n')
%constraints_W;






%NOT IMPORTANT
%{
for i = 1:m
    for j = 1:n
        if Tbin(i,j) == 0     
            %Uij = ;     % Performs the symbolic operation U_tilda*M_tilda+M*w GY(i,j)=\sum_l G(i,l)Y(l,j)
            [num,~] = numden(Y(i,j));
            cc      = coeffs(num,s);                                  % All elements of this vector must be 0....
            A_eq    = equationsToMatrix(cc,[vec(CQYs);vec(DQYs)]);      % Express system of equations in matrix form in terms of the vectorized versions of CQs and DQs
            A_eqs   = double(A_eq);    %A_eqs is the same as A_eq, for computation with sdpvars
            
            %for(f=1:size(A_equations,1)) %ALTERNATIVE FORM
            %    for(g=1:n*N)
            %Constraints=[Constraints, A_equations_normal(f,[(g-1)*m+1:g*m])*CQv(:,g)==0];
            %    end
            %end
            
            Constraints = [Constraints, A_eqs*[vec(CQYv);vec(DQYv)] == 0]; % Add the constraints in terms of the sdpvars CQv and DQv, by using A_eqs computed with symbolics
        end
    end
end

eq=eye(n)+Gs*Y;
for i = 1:m
    for j = 1:n
        if Rbin(i,j) == 0     
            %Uij = ;     % Performs the symbolic operation U_tilda*M_tilda+M*w GY(i,j)=\sum_l G(i,l)Y(l,j)
            [num,~] = numden(eq(i,j));
            cc      = coeffs(num,s);                                  % All elements of this vector must be 0....
            A_eq    = equationsToMatrix(cc,[vec(CQYs);vec(DQYs)]);      % Express system of equations in matrix form in terms of the vectorized versions of CQs and DQs
            A_eqs   = double(A_eq);    %A_eqs is the same as A_eq, for computation with sdpvars
            
            %for(f=1:size(A_equations,1)) %ALTERNATIVE FORM
            %    for(g=1:n*N)
            %Constraints=[Constraints, A_equations_normal(f,[(g-1)*m+1:g*m])*CQv(:,g)==0];
            %    end
            %end
            
            Constraints = [Constraints, A_eqs*[vec(CQYv);vec(DQYv)] == 0]; % Add the constraints in terms of the sdpvars CQv and DQv, by using A_eqs computed with symbolics
        end
    end
end


%}


%{
fprintf('Step 1-b: Encoding the constraint W=XG')
equation=W-X*Gs; %will set this to 0
for i = 1:n
    for j = 1:m
        fprintf('   Percentage %6.4f \n', 100*(m*(i-1)+j)/m/n );
        [num,~] = numden(equation(i,j));
        cc      = coeffs(num,s);                                  % All elements of this vector must be 0....
        [A_eq,b_eq]    = equationsToMatrix(cc,[vec(CQWs);vec(DQWs);vec(CQXs);vec(DQXs)]);      % Express system of equations in matrix form in terms of the vectorized versions of CQs and DQs
        A_eqs   = double(A_eq);    %A_eqs is the same as A_eq, for computation with sdpvars
        b_eqs    = double(b_eq);
        Constraints = [Constraints, A_eqs*[vec(CQWv);vec(DQWv);vec(CQXv);vec(DQXv)]== b_eqs]; % Add the constraints in terms of the sdpvars CQv and DQv, by using A_eqs computed with symbolics
    end
end
fprintf('Done \n')
fprintf('Step 1-c: Encoding the constraint Z=I+YG')
equation=Z-eye(m)-Y*Gs; %will set this to 0
for i = 1:m
    for j = 1:m
        fprintf('   Percentage %6.4f \n', 100*(m*(i-1)+j)/m/m );
        [num,~] = numden(equation(i,j));
        cc      = coeffs(num,s);                                  % All elements of this vector must be 0....
        [A_eq,b_eq]    = equationsToMatrix(cc,[vec(CQZs);vec(DQZs);vec(CQYs);vec(DQYs)]);      % Express system of equations in matrix form in terms of the vectorized versions of CQs and DQs
        A_eqs   = double(A_eq);    %A_eqs is the same as A_eq, for computation with sdpvars
        b_eqs    = double(b_eq);
        Constraints = [Constraints, A_eqs*[vec(CQZv);vec(DQZv);vec(CQYv);vec(DQYv)]== b_eqs]; % Add the constraints in terms of the sdpvars CQv and DQv, by using A_eqs computed with symbolics
    end
end

%{
%%%Sparsity constraints X \in Sparse(R), Y \in Sparse(T)
fprintf('Step 2-a: Encoding the constraint Y in Sparse(T) ...')
for i = 1:m
    for j = 1:n
        fprintf('   Percentage %6.4f \n', 100*(n*(i-1)+j)/m/n );
        
        if Tbin(i,j) == 0     
            [num,~] = numden(Y(i,j));
            cc      = coeffs(num,s);                                  % All elements of this vector must be 0....
            A_eq    = equationsToMatrix(cc,[vec(CQYs);vec(DQYs)]);      % Express system of equations in matrix form in terms of the vectorized versions of CQs and DQs
            A_eqs   = double(A_eq);    %A_eqs is the same as A_eq, for computation with sdpvars           
            Constraints = [Constraints, A_eqs*[vec(CQYv);vec(DQYv)] == 0]; % Add the constraints in terms of the sdpvars CQv and DQv, by using A_eqs computed with symbolics
        end
    end
end
fprintf('Done \n')
fprintf('Step 2-b: Encoding the constraint X in Sparse(R) ... only if not QI')


if QI == 0        %This cycle is useless if QI (redundant constraints). Hence, we skip it in this case.
    for i = 1:n
        for j = 1:n
            fprintf('   Percentage %6.4f \n', 100*(n*(i-1)+j)/n/n );
            if Rbin(i,j) == 0     % Whenever we need GY(i,j) = 0 ....
                [num,~] = numden(X(i,j));
                cc      = coeffs(num,s);                                  % All elements of this vector must be 0....
                A_eq    = equationsToMatrix(cc,[vec(CQXs);vec(DQXs)]);      % Express system of equations in matrix form in terms of the vectorized versions of CQs and DQs
                A_eqs   = double(A_eq);    %A_eqs is the same as A_eq, for computation with sdpvars        
                Constraints = [Constraints, A_eqs*[vec(CQv);vec(DQv)] == 0]; % Add the constraints in terms of the sdpvars CQv and DQv, by using A_eqs computed with symbolics
            end
        end
    end
end
fprintf('Done \n')
%}
%}
