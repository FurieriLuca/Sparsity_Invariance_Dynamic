clear all;
clc;


N=8;                                                                        % order of the controller
a=2;                                                                        % defines the basis for RH_infinity as {1/(s+a)^i}

generate_plant_data;                                                        % plant definition and relevant data
load('plant_data.mat');



%%QI Sparsity Constraints defined in [1]
Sbin1=[0 0 0 0 0;0 1 0 0 0;0 1 0 0 0;0 1 0 0 0;0 1 0 0 1];
Sbin2=[0 0 0 0 0;0 1 0 0 0;0 1 0 0 0;0 1 0 0 0;1 1 0 0 1];
Sbin3=[0 0 0 0 0;0 1 0 0 0;0 1 0 0 0;1 1 0 0 0;1 1 0 0 1];
Sbin4=[0 0 0 0 0;0 1 0 0 0;0 1 0 0 0;1 1 0 0 0;1 1 1 0 1];
Sbin5=[0 0 0 0 0;0 1 0 0 0;0 1 0 0 0;1 1 1 0 0;1 1 1 0 1];
Sbin6=[1 0 0 0 0;1 1 0 0 0;1 1 1 0 0;1 1 1 1 0;1 1 1 1 1];
centralized=ones(m,n);

%%non-QI Sparsity Constraint we study
Sbin7=[1 0 0 0 0;
       0 1 0 0 0;
       0 1 0 0 0;
       0 1 0 0 0;
       0 1 0 0 1]; 
Sbin=Sbin7;
QI=test_QI(Sbin,Delta)                                                      %variable QI used to avoid adding useless constraints later if QI=1 and reduce execution time

Tbin=Sbin;                                                                  % Matrix "T"                                                                 
Rbin=generate_SXlessS(Tbin);                                                % Matrix "R_MSI"
%Rbin=eye(n);

sparsity_constraints;                                                       % Creates matrices to encode sparsity constraints Y(s) \in Sparse(T) and G(s)Y(s) in Sparse(R) 

%%H2 norm minimization, SDP formulation
Qstatic=[CQv DQv];
P=sdpvar(size(A1_hat,1),size(A1_hat,1));
S=sdpvar(size(A1_hat,2),size(B2_hat,1));
R=sdpvar(size(A2_hat,1),size(A2_hat,1));
L=sdpvar(size(C2_hat,1),size(C2_hat,1));
gamma=sdpvar(1,1);

Constraints=[Constraints,trace(L)<=gamma, P>=0,R>=0];                                                                % (27)-(28) in our paper, Section V
Constraints=[Constraints, 
    [A1_hat*P+P*A1_hat' A1_hat*S-S*A2_hat+A_hat+B_hat*Qstatic*C_hat B1_hat+B_hat*Qstatic*F_hat-S*B2_hat;
    (A1_hat*S-S*A2_hat+A_hat+B_hat*Qstatic*C_hat)' R*A2_hat+A2_hat'*R R*B2_hat;
    (B1_hat+B_hat*Qstatic*F_hat-S*B2_hat)' B2_hat'*R -eye(size(B2_hat,2))]<=0];                                      % (29) 

Constraints=[Constraints, 
    [P zeros(size(P,1),size(R,2)) P*C1_hat';zeros(size(R,1),size(P,2)) R (C2_hat+E_hat*Qstatic*C_hat+C1_hat*S)';
    C1_hat*P C2_hat+E_hat*Qstatic*C_hat+C1_hat*S L]>=0,DQv==0];                                                      % (30), DQ=0 to guarantee that \mathcal{D}=0.

options=sdpsettings('allownonconvex',0,'solver','mosek','verbose',1);
sol=optimize(Constraints,gamma,options)

CQ=round(value(CQv),6);  %rounding to avoid false non-zeros                                                                                                                       
DQ=round(value(DQv),6);

%%RECOVER Q(s) and K(s) to check sparsities
for(i=1:n)
    for(j=1:m)
        Y(i,j)=CQ(i,[(j-1)*N+1:j*N])*inv(s*eye(size(AiQ,1))-AiQ)*BiQ+DQ(i,j);
    end
end
K=Y*inv(eye(n)+G*Y);

GQ=G*Y;
GQsubs=subs(GQ,s,rand);                                                     %just get rid of s to get the sparsity
GQbin=bin(GQsubs);  
Ksubs=subs(K,s,rand);                                                       %just get rid of s to get the sparsity
Kbin=bin(Ksubs);

sqrt(value(gamma))                                                          %value of the H2 norm!


%REFERENCES
%[1] ?"A Characterization of Convex Problems in Decentralized Control",
%[2]  "?Q-Parametrization and an SDP for H ?-optimal Decentralized Control"
%[3]  "?An efficient solution to multi-objective control problems with LMI
%objectives"