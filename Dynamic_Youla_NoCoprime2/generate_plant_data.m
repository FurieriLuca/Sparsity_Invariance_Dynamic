%PLANT DEFINITION
% Example taken from [1]

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

 % minimal state space realizations for P11, P12, P21
SYS = ss(P11,'min');      
A1  = SYS.A; B1  = SYS.B; C1  = SYS.C; D1  = SYS.D;

SYS = ss(P12,'min');
A2  = SYS.A; B2  = SYS.B; C2  = SYS.C; D2  = SYS.D;

SYS = ss(P21,'min');
A3  = SYS.A; B3  = SYS.B; C3  = SYS.C; D3  = SYS.D;


%% Parameterization of Y(s)
% Decomposition matrices onto basis, so that Cij (sI-AiQ)^-1 BiQ + Dij = Y(i,j)
shift = 1;
AiQ   = -a*eye(order) + diag(ones(order-abs(shift),1),shift);     
AiQx = -ax*eye(order) + diag(ones(order-abs(shift),1),shift);
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