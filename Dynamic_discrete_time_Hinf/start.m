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
clc;
%clc;

N = 3;           % order of the controller
a = 0;           % defines the basis for RH_infinity as {1/(s+a)^i}

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
Sbin6 = [1 0 0 0 0;1 1 0 0 0;1 1 1 0 0;1 1 1 1 0;1 1 1 1 1];

Sbin11 = [0 0 0 0 0;0 1 0 0 0;0 1 1 0 0;0 1 1 1 0;0 1 1 1 1]; %QI, closest to Sbin10, not feasible!

centralized = ones(m,n);

%% non-QI Sparsity Constraint we study
Sbin7 = [1 0 0 0 0;
    1 1 0 0 0;
    0 1 0 0 0;
    0 1 0 0 0;
    0 1 0 0 1];
Sbin = ones(5,5);
QI   = test_QI(Sbin,Delta);        % variable QI used to avoid adding useless constraints later if QI=1 and reduce execution time

Tbin = Sbin;                       % Matrix "T"
Rbin = generate_SXlessS(Tbin);     % Matrix "R_MSI"
%Rbin=eye(n);

fprintf('==================================================\n')
fprintf('              Sparisty Invariance Approch         \n')
fprintf('==================================================\n')

%% Encode sparsity Y(s) \in Sparse(T) 
Constraints = []
sparsity_constraints;

%% H2 norm minimization, SDP formulation
fprintf('Step 4: Encoding the other LMI constraint ...')
Qstatic = [CQYv DQYv];
E       = sdpvar(size(A1_hat,1),size(A1_hat,1));
R       = sdpvar(size(A2_hat,1),size(A2_hat,1));
S      = sdpvar(size(A1_hat,2),size(A2_hat,1));
gamma   = sdpvar(1,1);


%SDP constraints%
A_tilda=[A1_hat*E A1_hat*S+A_hat+B_hat*Qstatic*C_hat-S*A2_hat;
              zeros(size(A2_hat,1),size(A1_hat,2)) R*A2_hat];
B_tilda=[B1_hat+B_hat*Qstatic*F_hat-S*B2_hat;
              R*B2_hat];
C_tilda=[C1_hat*E C2_hat-E_hat*Qstatic*C_hat+C1_hat*S];
D_tilda=[D_hat-E_hat*Qstatic*F_hat];
X_tilda=blkdiag(E,R);

ABCD=[A_tilda B_tilda;
            C_tilda D_tilda];
 X_tilda_gamma=[X_tilda zeros(size(X_tilda,1),size(B_tilda,2));
                          zeros(size(B_tilda,2),size(X_tilda,2)) gamma*eye(size(B_tilda,2))];

big_matrix_constraints=[X_tilda_gamma ABCD';
                                      ABCD X_tilda_gamma];

Constraints = [Constraints,  big_matrix_constraints>=0];                                         % (17) in [1]                                         % (30), DQ=0 to guarantee that \mathcal{D}=0.

fprintf('Done \n')

% options = sdpsettings('allownonconvex',0,'solver','mosek','verbose',1);
fprintf('Step 5: call SDP solver to obtain a solution ... \n')
fprintf('=====================Yr=============================\n')

options = sdpsettings('allownonconvex',0,'solver','mosek','verbose',1);
sol     = optimize(Constraints,gamma,options);

Vgamma =value(gamma);  % value of the H2 norm!
fprintf('\n Hinf norm of the closed loop system is %6.4f \n', Vgamma);




%disp('\n\n\n We will now perform stability checks\n')
%checks;
