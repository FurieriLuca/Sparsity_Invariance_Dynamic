%%%%%%%
%This file is equivalent to start, but it runs the optimzation for all
%values of a in the in the interval [0.1,10] with a step of 0.1
%It generates a corresponding file "results.txt" with the results. The best
%choice is for a=2.

clear all;
clc;
fid = fopen('results.txt','w');

%clc;
for(i=1:100)
N = 6;           % order of the controller
a = i/10;   %2 %2%2%2%2        % defines the basis for RH_infinity as {1/(s+a)^i}

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
Sbin8 = [0 1 0 0 0;0 1 0 0 0;0 1 1 0 0;0 1 1 1 0;0 1 1 1 1]; %QI-closest to Sbin9
Sbin9 = [1 1 0 0 0;1 1 0 0 0;0 1 1 0 0;0 1 1 1 0;0 1 1 1 1]; %not QI
Sbin10 = [1 0 0 0 0;1 1 0 0 0;0 1 1 0 0;0 1 1 1 0;0 1 1 1 1]; %not QI, feasible with sparsity invariance
Sbin11 = [0 0 0 0 0;0 1 0 0 0;0 1 1 0 0;0 1 1 1 0;0 1 1 1 1]; %QI, closest to Sbin10, not feasible!

Sbin = Sbin7;
QI   = test_QI(Sbin,Delta);        % variable QI used to avoid adding useless constraints later if QI=1 and reduce execution time

Tbin = Sbin;                       % Matrix "T"
Rbin = generate_SXlessS(Tbin);     % Matrix "R_MSI"
%Rbin=eye(n);

fprintf('==================================================\n')
fprintf('              Sparsity Invariance Approch         \n')
fprintf('==================================================\n')

%% Encode sparsity Y(s) \in Sparse(T) and G(s)Y(s) in Sparse(R)
sparsity_constraints2;

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

options = sdpsettings('allownonconvex',0,'solver','mosek','verbose',1);
sol     = optimize(Constraints,gamma,options);

%CQ   = round(value(CQv),6);  %rounding to avoid false non-zeros
%DQ   = round(value(DQv),6);
CQ  = value(CQv);
DQ  = value(DQv);

Vgamma = sqrt(value(gamma));  % value of the H2 norm!
fprintf('\n H2 norm of the closed loop system is %6.4f \n', Vgamma);
fprintf(fid,'a=%f, norm=%f\n',a,Vgamma);

end

%{
%% RECOVER Q(s) and K(s) to check sparsities
Gi = (s*eye(size(AiQ,1))-AiQ)\BiQ;
for i = 1:n
    for j = 1:m
        Y(i,j) = CQ(i,(j-1)*N+1:j*N)*Gi + DQ(i,j);
    end
end
YQ=(Knom+inv(eye(5)-Knom*G)*Y)*inv(eye(5)-G*Knom);

K=YQ*inv(eye(5)+G*YQ);

Ksubs   = double(subs(K,s,rand));            %just get rid of s to get the sparsity
Kbin    = bin(Ksubs);

%K=simplify_tf(K);
pretty_K;

%GQ      = G*Y;
%GQsubs  = double(subs(GQ,s,rand));           %just get rid of s to get the sparsityKbin
%GQbin   = bin(GQsubs);
Ksubs   = double(subs(K,s,rand));            %just get rid of s to get the sparsity
Kbin    = bin(Ksubs);

for(i=1:m)
    for(j=1:n)
        if(Kbin(i,j)==0)
            K(i,j)=0;
        end
    end
end

Gmodel=syms2tf(G);
Kmodel=syms2tf(K_pretty);
Kmodel2=syms2tf(K);

%stability check
%Closed_Loop=K*inv(eye(n)-Gs*K)*P21s;
%for(i=1:m)
%    for(j=1:n)
%        if(isstable(syms2tf(Closed_Loop(i,j)))==0)
%            disp('PROBLEM: Closed Loop is not stable')
%        end
%    end
%end
%{
Closed_Loop=Ktf*inv(eye(n)-G*Ktf);
for(i=1:m)
    for(k=1:n)
        if(isstable(Closed_Loop(i,j))==0)
            i
            j
        end
    end
end
%}




%%  Verfication -- transfer functions
% s  = tf('s');
% Y = CQ*((s*eye(size(AiQ,1)*n)-kron(eye(n),AiQ))\kron(eye(n),BiQ)) + DQ;
% K = Y/(eye(n)+G*Y);
%
% Gdz = P11 + P12*(K/(eye(n) - G*K))*P21
% norm(Gdz,2)

%}

