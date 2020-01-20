clear all;
clc;

s=tf('s');
G=1/(s-1);  %plant


[Xl,Yl,Ur,Yr,Vl,Ul,Vr,Xr] = doubly_coprime_factorization(G);  %can modify poles of corresponding controller-observer


s=sym('s');


%%this is to verify that there are no unstable poles in our parameters, for
%%the chosen Y
fprintf('Consider the unstable plant \n')
Gs=1/(s-1)

s=tf('s');
G=1/(s-1);
SYS=ss(G,'min')
s=sym('s');

fprintf('We choose the paramater Y as\n')
Y=-4/(s+1)+8/(s+1)^2 %%computed so that X=1+G*Y and W=X*G are stable transfer functions

disp('which is such that X=1+GY, W=XG=G+GYG and Z=I+YG are stable:\n')
X=simplify(1+Gs*Y)
W=simplify(X*Gs)
Z=simplify(1+Y*Gs)

fprintf('We now want to verify that the corresponding Youla parameter is stable.\n The chosen doubly coprime factorization is\n')



Urs=(s-1)/(s+2)   %%Copied by hand by Ur, Yr, Ul computed above, in order to have symbolic variables. To try a different doubly-coprime parametrization, must modify values in the function doubly_coprime_factorization and copy results by hand
Yrs=-45/(s+2)
Uls=(s-1)/(s+14)

fprintf('resulting in the Youla parameter: \n')
Qs=simplify(Urs^-1*(Yrs-Y*Uls^-1)) %the resulting Youla parameter is stable. This is absolutely nontrivial from the expression for Qs. Now, we need to prove this...
