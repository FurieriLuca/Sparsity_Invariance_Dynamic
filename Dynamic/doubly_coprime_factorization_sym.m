function [M,N,M_tilda,N_tilda,U,V,U_tilda,V_tilda] = doubly_coprime_factorization_sym(G)
% Computes a doubly coprime factorization of G using "?A Connection Between
% State-Space and Doubly Coprime Fractional Representations" by C. N. Nett,
% C. A. Jacobson and M. J. Barls, TAC 1984

PLANT  = ss(G,'min');
A  = PLANT.A; B  = PLANT.B; C  = PLANT.C; D  = PLANT.D;
n  = size(A,1);
K  = place(A,B,-10000*ones(n,1));
F  = place(A',C',-10000*ones(n,1))';

s = sym('s');
N  = C*inv(s*eye(n)-A+B*K)*B;
M  = eye(n)-K*inv(s*eye(n)-A+B*K)*B;
U  = K*inv(s*eye(n)-A+F*C)*F;
V  = eye(n)+K*inv(s*eye(n)-A+F*C)*B;
M_tilda  = eye(n)-C*inv(s*eye(n)-A+F*C)*F;
N_tilda  = C*inv(s*eye(n)-A+F*C)*B;
V_tilda  = eye(n)+C*inv(s*eye(n)-A+B*K)*F;
U_tilda  = K*inv(s*eye(n)-A+B*K)*F;

end

