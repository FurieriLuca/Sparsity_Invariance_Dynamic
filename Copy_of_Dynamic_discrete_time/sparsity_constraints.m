



CXs = sym('CX',[n n*(N+1)]);
CXv = sdpvar(n,n*(N+1));          % decision variables

CYs = sym('CY',[m n*(N+1)]);
CYv = sdpvar(m,n*(N+1));          % decision variables


CWs = sym('CW',[n m*(N+1)]);
CWv = sdpvar(n,m*(N+1));          % decision variables


CZs = sym('CZ',[m m*(N+1)]);
CZv = sdpvar(m,m*(N+1));          % decision variables


%%Express X,Y,W,Z as FIR transfer matrices of order N
%{
shift=-1;
AiQ   =  diag(ones(N-abs(shift),1),shift);
AQ    = kron(eye(n),AiQ); %from size of P21...2n
I     = eye(N);
BiQ   = I(:,1);
BQ    = kron(eye(n),BiQ);
%}
X=zeros(n,n);
Y=zeros(m,n);
W=zeros(n,m);
Z=zeros(m,m);
for(t=1:N+1)
    X=X+CXs(:,[(t-1)*n+1:t*n])/z^(t-1);
    Y=Y+CYs(:,[(t-1)*n+1:t*n])/z^(t-1);
    W=W+CWs(:,[(t-1)*m+1:t*m])/z^(t-1);
    Z=Z+CZs(:,[(t-1)*m+1:t*m])/z^(t-1);
end

%{
for(i=1:n)
    for(j=1:n)
        X(i,j)=CXs(i,[(j-1)*N+1:j*N])*inv(z*eye(N)-AiQ)*BiQ+DXs(i,j);
    end
end
for(i=1:m)
    for(j=1:n)
        Y(i,j)=CYs(i,[(j-1)*N+1:j*N])*inv(z*eye(N)-AiQ)*BiQ+DYs(i,j);
    end
end
for(i=1:n)
    for(j=1:m)
        W(i,j)=CWs(i,[(j-1)*N+1:j*N])*inv(z*eye(N)-AiQ)*BiQ+DWs(i,j);
    end
end
for(i=1:m)
    for(j=1:m)
        Z(i,j)=CZs(i,[(j-1)*N+1:j*N])*inv(z*eye(N)-AiQ)*BiQ+DZs(i,j);
    end
end
%}
%%%achievability constraints
ach1 = X-G*Y-eye(n);
ach2 = W-G*Z;
ach3 = -X*G+W;
ach4 = -Y*G+Z-eye(m);
Constraints=[];


for(i=1:n)
    for(j=1:n)
        fprintf('   Percentage %6.4f \n', 100*(n*(i-1)+j)/n/n );
        [num,~]=numden(ach1(i,j));
        cc=coeffs(num,z);
        [A_eq,b_eq]    = equationsToMatrix(cc,[vec(CXs);vec(CYs)]);
        A_eqs   = double(A_eq);
        b_eqs    = double(b_eq);
        Constraints = [Constraints, A_eqs*[vec(CXv);vec(CYv)]== b_eqs];
    end
end

for(i=1:n)
    for(j=1:m)
        fprintf('   Percentage %6.4f \n', 100*(m*(i-1)+j)/m/n );
        [num,~]=numden(ach2(i,j));
        cc=coeffs(num,z);
        [A_eq,b_eq]    = equationsToMatrix(cc,[vec(CWs);vec(CZs)]);
        A_eqs   = double(A_eq);
        b_eqs    = double(b_eq);
        Constraints = [Constraints, A_eqs*[vec(CWv);vec(CZv)]== b_eqs];
    end
end
for(i=1:n)
    for(j=1:m)
        fprintf('   Percentage %6.4f \n', 100*(m*(i-1)+j)/m/n );
        [num,~]=numden(ach3(i,j));
        cc=coeffs(num,z);
        [A_eq,b_eq]    = equationsToMatrix(cc,[vec(CXs);vec(CWs)]);
        A_eqs   = double(A_eq);
        b_eqs    = double(b_eq);
        Constraints = [Constraints, A_eqs*[vec(CXv);vec(CWv)]== b_eqs];
    end
end
for(i=1:m)
    for(j=1:m)
        fprintf('   Percentage %6.4f \n', 100*(m*(i-1)+j)/m/m );
        [num,~]=numden(ach4(i,j));
        cc=coeffs(num,z);
        [A_eq,b_eq]    = equationsToMatrix(cc,[vec(CYs);vec(CZs)]);
        A_eqs   = double(A_eq);
        b_eqs    = double(b_eq);
        Constraints = [Constraints, A_eqs*[vec(CYv);vec(CZv)]== b_eqs];
    end
end

%% Constraint "Achievability": X=I+GY,W=XG,Z=I+YG
%{
fprintf('Step 1-a: Encoding the constraint X=I+GY is stable \n')
constraints_X;


fprintf('Step 1-b: Encoding the constraint W=(I+GY)G is stable \n')
constraints_W;
%}





%NOT IMPORTANT
%{
for i = 1:m
    for j = 1:n
        if Tbin(i,j) == 0
            %Uij = ;     % Performs the symbolic operation U_tilda*M_tilda+M*w GY(i,j)=\sum_l G(i,l)Y(l,j)
            [num,~] = numden(Y(i,j));
            cc      = coeffs(num,s);                                  % All elements of this vector must be 0....
            A_eq    = equationsToMatrix(cc,[vec(CYs);vec(DYs)]);      % Express system of equations in matrix form in terms of the vectorized versions of CQs and DQs
            A_eqs   = double(A_eq);    %A_eqs is the same as A_eq, for computation with sdpvars
            
            %for(f=1:size(A_equations,1)) %ALTERNATIVE FORM
            %    for(g=1:n*N)
            %Constraints=[Constraints, A_equations_normal(f,[(g-1)*m+1:g*m])*CQv(:,g)==0];
            %    end
            %end
            
            Constraints = [Constraints, A_eqs*[vec(CYv);vec(DYv)] == 0]; % Add the constraints in terms of the sdpvars CQv and DQv, by using A_eqs computed with symbolics
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
            A_eq    = equationsToMatrix(cc,[vec(CYs);vec(DYs)]);      % Express system of equations in matrix form in terms of the vectorized versions of CQs and DQs
            A_eqs   = double(A_eq);    %A_eqs is the same as A_eq, for computation with sdpvars
            
            %for(f=1:size(A_equations,1)) %ALTERNATIVE FORM
            %    for(g=1:n*N)
            %Constraints=[Constraints, A_equations_normal(f,[(g-1)*m+1:g*m])*CQv(:,g)==0];
            %    end
            %end
            
            Constraints = [Constraints, A_eqs*[vec(CYv);vec(DYv)] == 0]; % Add the constraints in terms of the sdpvars CQv and DQv, by using A_eqs computed with symbolics
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
        [A_eq,b_eq]    = equationsToMatrix(cc,[vec(CWs);vec(DWs);vec(CXs);vec(DXs)]);      % Express system of equations in matrix form in terms of the vectorized versions of CQs and DQs
        A_eqs   = double(A_eq);    %A_eqs is the same as A_eq, for computation with sdpvars
        b_eqs    = double(b_eq);
        Constraints = [Constraints, A_eqs*[vec(CWv);vec(DWv);vec(CXv);vec(DXv)]== b_eqs]; % Add the constraints in terms of the sdpvars CQv and DQv, by using A_eqs computed with symbolics
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
        [A_eq,b_eq]    = equationsToMatrix(cc,[vec(CZs);vec(DZs);vec(CYs);vec(DYs)]);      % Express system of equations in matrix form in terms of the vectorized versions of CQs and DQs
        A_eqs   = double(A_eq);    %A_eqs is the same as A_eq, for computation with sdpvars
        b_eqs    = double(b_eq);
        Constraints = [Constraints, A_eqs*[vec(CZv);vec(DZv);vec(CYv);vec(DYv)]== b_eqs]; % Add the constraints in terms of the sdpvars CQv and DQv, by using A_eqs computed with symbolics
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
            A_eq    = equationsToMatrix(cc,[vec(CYs);vec(DYs)]);      % Express system of equations in matrix form in terms of the vectorized versions of CQs and DQs
            A_eqs   = double(A_eq);    %A_eqs is the same as A_eq, for computation with sdpvars
            Constraints = [Constraints, A_eqs*[vec(CYv);vec(DYv)] == 0]; % Add the constraints in terms of the sdpvars CQv and DQv, by using A_eqs computed with symbolics
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
                A_eq    = equationsToMatrix(cc,[vec(CXs);vec(DXs)]);      % Express system of equations in matrix form in terms of the vectorized versions of CQs and DQs
                A_eqs   = double(A_eq);    %A_eqs is the same as A_eq, for computation with sdpvars
                Constraints = [Constraints, A_eqs*[vec(CQv);vec(DQv)] == 0]; % Add the constraints in terms of the sdpvars CQv and DQv, by using A_eqs computed with symbolics
            end
        end
    end
end
fprintf('Done \n')
%}
%}
