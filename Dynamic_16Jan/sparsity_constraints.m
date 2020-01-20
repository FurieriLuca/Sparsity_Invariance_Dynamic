
s = sym('s');                                        % Redefines "s" as symbolic variable instead of a tf variable, which would be wrongly interpreted by matlab later
assume(s,'real');
Gs = 1/(s-0.1)  
Knoms = [0 0 0 0 0;0 -6/(s+3) 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 0 0 -6/(s+3)];       % nominal stabilizing controller in the structure
%asd=load('K2.mat');
%Knoms=asd.K2;
%Gs=Gs*inv(eye(m)-Knoms*Gs);




% Defines CQ and DQ both as symbolic variables and sdpvar variables, for
% four paramenters X=I+GY,Y,W=(I+GY)G,Z=I+YG

Gs;

CQYs = sym('CQY',[m n*N]);
DQYs = sym('DQY',[m n]);
CQYv = sdpvar(m,n*N);          % decision variables
DQYv = sdpvar(m,n);



%Force Y in RHinfinity

for(i=1:m)
    for(j=1:n)
        Y(i,j)=(CQYs(i,[(j-1)*N+1:j*N])*inv(s*eye(N)-AiQ)*BiQ+DQYs(i,j));
    end
end



Constraints = [];

%% Constraint "Achievability": X=I+GY,W=XG,Z=I+YG

fprintf('Step 1-a: Encoding the constraint X=I+GY is stable \n')
%constraints_X;
  


fprintf('Step 1-b: Encoding the constraint W=(I+GY)G is stable \n')
constraints_W;


%{
fprintf('Step 1-c: Encoding the constraint X=I+GY is stable \n')
constraints_Z;
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
