function [Xl,Yl,Ur,Yr,Vl,Ul,Vr,Xr] = doubly_coprime_factorization(G)
% Computes a doubly coprime factorization of G using "?A Connection Between
% State-Space and Doubly Coprime Fractional Representations" by C. N. Nett,
% C. A. Jacobson and M. J. Barls, TAC 1984

s = tf('s');
PLANT  = ss(G,'min');
A  = PLANT.A; B  = PLANT.B; C  = PLANT.C; D  = PLANT.D;
n  = size(A,1);
K  = place(A,B,-2*ones(n,1));
F  = place(A',C',-14*ones(n,1))';

N  = C*inv(s*eye(n)-A+B*K)*B;
M  = eye(n)-K*inv(s*eye(n)-A+B*K)*B;
U  = K*inv(s*eye(n)-A+F*C)*F;
V  = eye(n)+K*inv(s*eye(n)-A+F*C)*B;
M_tilda  = eye(n)-C*inv(s*eye(n)-A+F*C)*F;
N_tilda  = C*inv(s*eye(n)-A+F*C)*B;
V_tilda  = eye(n)+C*inv(s*eye(n)-A+B*K)*F;
U_tilda  = K*inv(s*eye(n)-A+B*K)*F;

Xl=V;
Yl=-U;
Ur=M;
Yr=-U_tilda;
Vl=N_tilda;
Ul=M_tilda;
Vr=N;
Xr=V_tilda;

end

