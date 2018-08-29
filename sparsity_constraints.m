s=sym('s');                                                                                 %Redefines "s" as symbolic variable instead of a tf variable, which would be wrongly interpreted by matlab later
G=[1/(s+1) 0 0 0 0;                                                                         %Redefines the plant in terms of symbolic "s"
    1/(s+1) 1/(s-1) 0 0 0;
    1/(s+1) 1/(s-1) 1/(s+1) 0 0;
    1/(s+1) 1/(s-1) 1/(s+1) 1/(s+1) 0;
    1/(s+1) 1/(s-1) 1/(s+1) 1/(s+1) 1/(s-1)];

I=eye(n);

CQs=sym('CQ',[m n*N]);                                                                      %Defines CQ and DQ both as symbolic variables and sdpvar variables
DQs=sym('DQ',[m n]);
CQv=sdpvar(m,n*N);
DQv=sdpvar(m,n);

Constraints=[];


%Creates the constraint Y(s) \in Sparse(T) in terms of CQv and DQv (same as [2])
for(i=1:m)                                                                                  
    for(j=1:n)
        if(Tbin(i,j)==0) %main cycle
            Constraints=[Constraints, CQv(i,[(j-1)*N+1:j*N])==0, DQv(i,j)==0];
        end
    end
end

%Creates the constraint G(s)Y(s) \in Sparse(R) in terms of CQv and DQv
if(QI==0)                                                                                   %This cycle is useless if QI (redundant constraints). Hence, we skip it in this case.
    for(i=1:n)                                                                                                                
        for(j=1:n)
            completion=100*(n*(i-1)+j)/n/n                                                  % Just shows progress in the interval 0-100
            
            if(Rbin(i,j)==0 )                                                               % Whenever we need GY(i,j)=0....
                sum=0;
                for(l=1:m) 
                    Qslj=CQs(l,[(j-1)*N+1:j*N])*inv(s*eye(N)-AiQ)*BiQ+DQs(l,j);
                    sum=sum+G(i,l)*Qslj;                                                    %Performs the symbolic matrix product GY(i,j)=\sum_l G(i,l)Y(l,j)
                end
                [num,den]=numden(sum);
                cc=coeffs(num,s);                                                           %All elements of this vector must be 0....
                A_eq=equationsToMatrix(cc,[vec(CQs);vec(DQs)]);                             %Express system of equations in matrix form in terms of the vectorized versions of CQs and DQs
                
                for(f=1:size(A_eq,1))                                                       %A_eqs is the same as A_eq, just rewritten as a standard matrix which allows use with sdpvars
                    for(g=1:size(A_eq,2))
                        A_eqs(f,g)=double(A_eq(f,g));
                    end
                end
                
%for(f=1:size(A_equations,1)) %ALTERNATIVE FORM
%    for(g=1:n*N)
%Constraints=[Constraints, A_equations_normal(f,[(g-1)*m+1:g*m])*CQv(:,g)==0];
%    end
%end
                Constraints=[Constraints, A_eqs*[vec(CQv);vec(DQv)]==0];                    % Add the constraints in terms of the sdpvars CQv and DQv, by using A_eqs computed with symbolics
            end
        end
    end
end
