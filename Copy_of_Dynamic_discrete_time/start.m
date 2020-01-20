% Solve an optimal decentralized control problem using the sparisty invariance approch [0]
%
%         min <K> ||P11 + P12 * K(I - GK)^(-1) * P21||
%         s.t.      K \in S
%
% where P11, P12, G, P21 are transfer functions (plant dynamics), and S is
% a given sparsity constraint. In [0], we propose to retrict the problem into
%
%         min <X,Y> ||P11 + P12 * Y * P21||
%         s.t.      Y, X, XG, I+YG \in RH
%                   X - I = GY
%                   Y \in T, X \in R
%
% where RH denotes real-rational proper stable transfer functions,
% T, R are sparsity constraints, coming from the sparsity invirance approach
%
% Authors: L. Furieri, Automatic Control Laboratory, ETH.
%          Y. Zheng, Department of Engineering Science, University of Oxford
%
% References
% [0] "Convex Restrictions in Distributed Control"
% [1] "A Characterization of Convex Problems in Decentralized Control",
% [2] "Q-Parametrization and an SDP for H ?-optimal Decentralized Control"
% [3] "An efficient solution to multi-objective control problems with LMI objectives"

clear all;
clc;
%clc;

N = 7;


%% generate plant data
generate_plant_data;         % plant definition and relevant data
%load('plant_data.mat');

%% QI Sparsity Constraints defined in [1]
Sbin1 = [0 0 0 0 0;0 1 0 0 0;0 1 0 0 0;0 1 0 0 0;0 1 0 0 1];
Sbin2 = [0 0 0 0 0;0 1 0 0 0;0 1 0 0 0;0 1 0 0 0;1 1 0 0 1];
Sbin3 = [0 0 0 0 0;0 1 0 0 0;0 1 0 0 0;1 1 0 0 0;1 1 0 0 1];
Sbin4 = [0 0 0 0 0;0 1 0 0 0;0 1 0 0 0;1 1 0 0 0;1 1 1 0 1];
Sbin5 = [0 0 0 0 0;0 1 0 0 0;0 1 0 0 0;1 1 1 0 0;1 1 1 0 1];
Sbin6 = [1 0 0 0 0;1 1 0 0 0;1 1 1 0 0;1 1 1 1 0;1 1 1 1 1];
centralized = ones(m,n);

%% non-QI Sparsity Constraint we study
Sbin7 = [1 0 0 0 0;
    1 1 0 0 0;
    0 1 0 0 0;
    0 1 0 0 0;
    0 1 0 0 1];
Sbin = Sbin1;
QI   = test_QI(Sbin,Delta);        % variable QI used to avoid adding useless constraints later if QI=1 and reduce execution time

Tbin = Sbin;                       % Matrix "T"
Rbin = generate_SXlessS(Tbin);     % Matrix "R_MSI"
%Rbin=eye(n);

fprintf('==================================================\n')
fprintf('              Sparsity Invariance Approch         \n')
fprintf('==================================================\n')

%% Encode sparsity Y(s) \in Sparse(T) and G(s)Y(s) in Sparse(R)
sparsity_constraints;

%% Hinf norm minimization, SDP formulation
fprintf('Step 4: Encoding the other LMI constraint ...')

%{
Qstatic = [CQYv DQYv];
S       = sdpvar(size(A1_hat,2),size(A2_hat,1));
R       = sdpvar(size(A2_hat,1),size(A2_hat,1));
E       = sdpvar(size(A1_hat,2),size(A1_hat,2));
gamma   = sdpvar(1,1);

%%EQ (18) in [2]
A_tilda = [A1_hat*E A1_hat*S+A_hat+B_hat*Qstatic*C_hat-S*A2_hat;
           zeros(size(R,1),size(E,2)) R*A2_hat];
B_tilda = [B1_hat+B_hat*Qstatic*F_hat-S*B2_hat;
           R*B2_hat];
C_tilda = [C1_hat*E C2_hat-E_hat*Qstatic*C_hat+C1_hat*S];
D_tilda = [D_hat-E_hat*Qstatic*F_hat];
X_tilda = blkdiag(E,R);

%additional block matrices for compactness
X_tilda_gamma = [X_tilda zeros(size(X_tilda,1),size(B_tilda,2));
               zeros(size(B_tilda,2),size(X_tilda,2)) gamma*eye(size(B_tilda,2))];
ABCD_tilda = [A_tilda B_tilda;
              C_tilda D_tilda];

epsilon = 1e-4;
I = eye(size(X_tilda,1)+size(B_tilda,2)+size(A_tilda,1)+size(C_tilda,1));

Constraints=[];
Constraints = [Constraints,
               [X_tilda_gamma ABCD_tilda';
                ABCD_tilda X_tilda_gamma] >= 0];                                        

fprintf('Done \n')
%}
% options = sdpsettings('allownonconvex',0,'solver','mosek','verbose',1);

%%H2cost
cost_matrix=[CWv(:,[(1-1)*m+1:1*m]) CXv(:,[(1-1)*n+1:1*n])-eye(n,n);CZv(:,[(1-1)*m+1:1*m])-eye(m) CYv(:,[(1-1)*n+1:1*n])]'*[CWv(:,[(1-1)*m+1:1*m]) CXv(:,[(1-1)*n+1:1*n])-eye(n,n);CZv(:,[(1-1)*m+1:1*m])-eye(m) CYv(:,[(1-1)*n+1:1*n])];
cost=trace(cost_matrix);
%%build(W[t])
for(t=2:N+1)
       cost_matrix=[CWv(:,[(t-1)*m+1:t*m]) CXv(:,[(t-1)*n+1:t*n])-eye(n,n);CZv(:,[(t-1)*m+1:t*m])-eye(m) CYv(:,[(t-1)*n+1:t*n])]'*[CWv(:,[(t-1)*m+1:t*m]) CXv(:,[(t-1)*n+1:t*n])-eye(n,n);CZv(:,[(t-1)*m+1:t*m])-eye(m) CYv(:,[(t-1)*n+1:t*n])];
       cost=cost+trace(cost_matrix);
end


fprintf('=====================Yr=============================\n')

options = sdpsettings('allownonconvex',0,'solver','mosek','verbose',1);
sol     = optimize(Constraints,cost,options);

Vgamma = value(cost);  % value of the H2 norm!
fprintf('\n H2 norm of the closed loop system is %6.4f \n', sqrt(Vgamma));



%CHECKS
z=sym('z');
Gz = [0.1/(z-0.5) 0 0 0 0;                             
     0.1/(z-0.5) 1/(z-2) 0 0 0;
     0.1/(z-0.5) 1/(z-2) 0.1/(z-0.5) 0 0;                                              % plant tf definition
     0.1/(z-0.5) 1/(z-2) 0.1/(z-0.5) 0.1/(z-0.5) 0;
     0.1/(z-0.5) 1/(z-2) 0.1/(z-0.5) 0.1/(z-0.5) 1/(z-2)]; 
Xr=zeros(n,n);
Yr=zeros(m,n);
Wr=zeros(n,m);
Zr=zeros(m,m);
for(t=1:N+1)
    Xr=Xr+value(CXv(:,[(t-1)*n+1:t*n]))/z^(t-1);
    Yr=Yr+value(CYv(:,[(t-1)*n+1:t*n]))/z^(t-1);
    Wr=Wr+value(CWv(:,[(t-1)*m+1:t*m]))/z^(t-1);
    Zr=Zr+value(CZv(:,[(t-1)*m+1:t*m]))/z^(t-1);
end




%{
z=tf('z');
CQYr  = value(CQYv);
DQYr  = value(DQYv);
Gi = (z*eye(size(AiQ,1))-AiQ)\BiQ;
for(i=1:m)
    for(j=1:n)
       % Yr(i,j) = CQYr(i,(j-1)*N+1:j*N)*Gi + DQYr(i,j);
        Yr(i,j) = (z-2)/z^2;
    end
end
w_to_z=[G+G*Yr*G G*Yr;Yr*G eye(size(Yr,1),size(Yr,2))];
%}
