%PLANT DEFINITION
% Example took from 
% "A Characterization of Convex Problems in Decentralized Control",

s = tf('s');
G = [1/(s+1) 0 0 0 0;                             
     1/(s+1) 1/(s-1) 0 0 0;
     1/(s+1) 1/(s-1) 1/(s+1) 0 0;                                              % plant tf definition
     1/(s+1) 1/(s-1) 1/(s+1) 1/(s+1) 0;
     1/(s+1) 1/(s-1) 1/(s+1) 1/(s+1) 1/(s-1)]; 
Delta = [1 0 0 0 0;1 1 0 0 0;1 1 1 0 0;1 1 1 1 0;1 1 1 1 1];                   % plant binary structure
n     = size(G,1);
m     = size(G,2);
P11   = [G zeros(n,n);zeros(n,n) zeros(n,n)];                                   
P12   = [G;eye(n)];
P21   = [G eye(n)];

Knom = [0 0 0 0 0;0 -6/(s+3) 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 0 0 -6/(s+3)];       % nominal stable and stabilizing controller in the structure


% From Theorem 17 in [1]

M_tilda  = inv(eye(n)-G*Knom);
M  = -inv(eye(m)-Knom*G);
N_tilda  = G*inv(eye(m)-Knom*G);
N  = -G*inv(eye(m)-Knom*G);
V  = -eye(m);
U  =  Knom;
V_tilda  = eye(n);
U_tilda  = -Knom;


%[M,N,M_tilda,N_tilda,U,V,U_tilda,V_tilda]=doubly_coprime_factorization(G);
% see equation (3) in [1]
%T3   = [M_tilda N_tilda];                                                     
%T2   = [N;M];                                   
%T1   = [V_tilda*M_tilda -V_tilda*N_tilda;U_tilda*M_tilda eye(n)-U_tilda*N_tilda];
E1   = P11-P12*U_tilda*M_tilda*P21;
E2   = P12*M;
E3   = M_tilda*P21;
%T3   = (eye(n)-G*Knom)\P21;                                                     
%T1   = P11 + P12*Knom*T3;                                   
%T2   = -P12/(eye(m)-Knom*G);


 % minimal state space realizations for T1, T2, T3
SYS = ss(E1,'min');      
A1  = SYS.A; B1  = SYS.B; C1  = SYS.C; D1  = SYS.D;
disp('sonqua1')

SYS = ss(E2,'min');
A2  = SYS.A; B2  = SYS.B; C2  = SYS.C; D2  = SYS.D;
disp('sonqua2')

SYS = ss(E3,'min');
A3  = SYS.A; B3  = SYS.B; C3  = SYS.C; D3  = SYS.D;
disp('sonqua3')


%% Parameterization of Y(s)
% Decomposition matrices onto basis, so that Cij (sI-AiQ)^-1 BiQ + Dij = Y(i,j)
shift = 1;
AiQ   = -a*eye(order) + diag(ones(order-abs(shift),1),shift);        
AQ    = kron(eye(n),AiQ);
I     = eye(order);
BiQ   = I(:,order);
BQ    = kron(eye(n),BiQ);

%% State-space realization of the augmented LTI system
%  See (19) in [2] for the following matrices
A1_hat = blkdiag(A1, A2);
A2_hat = [A3 zeros(size(A3,1),size(AQ,2));BQ*C3 AQ];
A_hat  = zeros(size(A1_hat,1),size(A2_hat,2));
B1_hat = [B1;zeros(size(B2,1),size(B1,2))];
B_hat  = [zeros(size(B1,1),size(B2,2));B2];
B2_hat = [B3;BQ*D3];
C1_hat = [C1 -C2];
C_hat  = [zeros(n*order,size(C3,2)) eye(n*order);C3 zeros(size(C3,1),n*order)];
C2_hat = [zeros(size(C1_hat,1),size(C_hat,2))];
D_hat  = D1;
E_hat  = -D2;
F_hat  = [zeros(size(C_hat,1)-size(D3,1),size(D_hat,2));D3];

save plant_data