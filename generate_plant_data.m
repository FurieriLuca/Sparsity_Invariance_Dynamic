%PLANT DEFINITION
s=tf('s');
G=[1/(s+1) 0 0 0 0;                             
    1/(s+1) 1/(s-1) 0 0 0;
    1/(s+1) 1/(s-1) 1/(s+1) 0 0;                                            %plant tf definition
    1/(s+1) 1/(s-1) 1/(s+1) 1/(s+1) 0;
    1/(s+1) 1/(s-1) 1/(s+1) 1/(s+1) 1/(s-1)];                               
Delta=[1 0 0 0 0;1 1 0 0 0;1 1 1 0 0;1 1 1 1 0;1 1 1 1 1];                  %plant binary structure
n=size(G,1);
m=size(G,2);
P11=[G zeros(n,n);zeros(n,n) zeros(n,n)];                                   
P12=[G;eye(n)];
P21=[G eye(n)];

Knom=[0 0 0 0 0;0 -6/(s+3) 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 0 0 -6/(s+3)];     %nominal stabilizing controller in the structure
T1=P11+P12*Knom*inv(eye(n)-G*Knom)*P21;                                     %see equation (3) in [1]
T2=-P12*inv(eye(n)-Knom*G);
T3=inv(eye(n)-G*Knom)*P21;

SYS=ss(T1,'min')                                                            %minimal state space realizations for T1, T2, T3
A1=SYS.A;
B1=SYS.B;
C1=SYS.C;
D1=SYS.D;

SYS=ss(T2,'min')
A2=SYS.A;
B2=SYS.B;
C2=SYS.C;
D2=SYS.D;

SYS=ss(T3,'min')
A3=SYS.A;
B3=SYS.B;
C3=SYS.C;
D3=SYS.D;


shift=1;
AiQ=-a*eye(N)+diag(ones(N-abs(shift),1),shift);                             %Decomposition matrices onto basis, so that Cij(sI-AiQ)^-1BiQ+Dij=Y(i,j)
AQ=kron(eye(n),AiQ);
I=eye(N);
BiQ=I(:,N);
BQ=kron(eye(n),BiQ);

A1_hat=[A1 zeros(size(A1,1),size(A2,2));zeros(size(A2,1),size(A1,2)) A2];   % See (19) in [2] for the following matrices
A2_hat=[A3 zeros(size(A3,1),size(AQ,2));BQ*C3 AQ];
A_hat=zeros(size(A1_hat,1),size(A2_hat,2));
B1_hat=[B1;zeros(size(B2,1),size(B1,2))];
B_hat=[zeros(size(B1,1),size(B2,2));B2];
B2_hat=[B3;BQ*D3];
C1_hat=[C1 -C2];
C_hat=[zeros(n*N,size(C3,2)) eye(n*N);C3 zeros(size(C3,1),n*N)];
C2_hat=[zeros(size(C1_hat,1),size(C_hat,2))];
D_hat=D1;
E_hat=-D2;
F_hat=[zeros(size(C_hat,1)-size(D3,1),size(D_hat,2));D3];

save plant_data