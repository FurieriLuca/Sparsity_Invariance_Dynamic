% Solve an optimal decentralized control problem using the sparisty invirance approch [0]
%
%         min <K> ||P11 + P12 * K(I - GK)^(-1) * P21||                          
%         s.t.      K \in S
%
% where P11, P12, G, P21 are transfer functions (plant dynamics), and S is 
% a given sparsity constraint. In [0], we propose to retrict the problem into
%
%         min <X,Y> ||T1 - T2*Y*T3|| 
%         s.t.      Y \in RH
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
%clc;

N = 7;           % order of the controller
a = 2;           % defines the basis for RH_infinity as {1/(s+a)^i}

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
Sbin = Sbin7;
QI   = test_QI(Sbin,Delta);        % variable QI used to avoid adding useless constraints later if QI=1 and reduce execution time

Tbin = Sbin;                       % Matrix "T"                                                                 
Rbin = generate_SXlessS(Tbin);     % Matrix "R_MSI"
%Rbin=eye(n);

fprintf('==================================================\n')
fprintf('              Sparisty Invariance Approch         \n')
fprintf('==================================================\n')

%% Encode sparsity Y(s) \in Sparse(T) and G(s)Y(s) in Sparse(R)
sparsity_constraints;              

%% H2 norm minimization, SDP formulation
fprintf('Step 3: Encoding the other LMI constraint ...')
Qstatic = [CQv DQv];
P       = sdpvar(size(A1_hat,1),size(A1_hat,1));
S       = sdpvar(size(A1_hat,2),size(B2_hat,1));
R       = sdpvar(size(A2_hat,1),size(A2_hat,1));
L       = sdpvar(size(C2_hat,1),size(C2_hat,1));
gamma   = sdpvar(1,1);

Constraints = [Constraints,trace(L)<=gamma, P>=0,R>=0];                                                                   % (27)-(28) in our paper, Section V
Constraints = [Constraints, 
    [A1_hat*P+P*A1_hat' A1_hat*S-S*A2_hat+A_hat+B_hat*Qstatic*C_hat B1_hat+B_hat*Qstatic*F_hat-S*B2_hat;
    (A1_hat*S-S*A2_hat+A_hat+B_hat*Qstatic*C_hat)' R*A2_hat+A2_hat'*R R*B2_hat;
    (B1_hat+B_hat*Qstatic*F_hat-S*B2_hat)' B2_hat'*R -eye(size(B2_hat,2))] <= 0];                                         % (29) 

Constraints = [Constraints, 
    [P zeros(size(P,1),size(R,2)) P*C1_hat';zeros(size(R,1),size(P,2)) R (C2_hat+E_hat*Qstatic*C_hat+C1_hat*S)';
    C1_hat*P C2_hat+E_hat*Qstatic*C_hat+C1_hat*S L] >= 0,DQv == 0];                                                       % (30), DQ=0 to guarantee that \mathcal{D}=0.
fprintf('Done \n')

% options = sdpsettings('allownonconvex',0,'solver','mosek','verbose',1);
fprintf('Step 4: call SDP solver to obtain a solution ... \n')
fprintf('==================================================\n')

options = sdpsettings('allownonconvex',0,'solver','sedumi','verbose',1);
sol     = optimize(Constraints,gamma,options);

%CQ   = round(value(CQv),6);  %rounding to avoid false non-zeros                                                                                                                       
%DQ   = round(value(DQv),6);
CQ  = value(CQv);
DQ  = value(DQv);

Vgamma = sqrt(value(gamma));  % value of the H2 norm!
fprintf('\n H2 norm of the closed loop system is %6.4f \n', Vgamma);

%% RECOVER Q(s) and K(s) to check sparsities
Gi = (s*eye(size(AiQ,1))-AiQ)\BiQ;
for i = 1:n
    for j = 1:m
        Y(i,j) = CQ(i,(j-1)*N+1:j*N)*Gi + DQ(i,j);
    end
end
K = Y/(eye(n)+Gs*Y);
GQ      = Gs*Y;
GQsubs  = double(subs(GQ,s,rand));           %just get rid of s to get the sparsityKbin
GQbin   = bin(GQsubs);  
Ksubs   = double(subs(K,s,rand));            %just get rid of s to get the sparsity
Kbin    = bin(Ksubs);

%%  Verfication -- transfer functions
% s  = tf('s');
% Y = CQ*((s*eye(size(AiQ,1)*n)-kron(eye(n),AiQ))\kron(eye(n),BiQ)) + DQ;
% K = Y/(eye(n)+G*Y);
% 
% Gdz = P11 + P12*(K/(eye(n) - G*K))*P21
% norm(Gdz,2)



