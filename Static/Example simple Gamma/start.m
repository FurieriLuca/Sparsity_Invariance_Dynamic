


%%% Finding a Gamma to improve performance. Still NOT COHERENT

clear;close all

%A=3*[2 1 0;0 -5 1;-0.67 0.89 0.5];
%=10*[2 1 1;2 -5 1;1 -0.89 0.5];
A=[2 1 5;0 -1 1;-1 1 0.5];
%B2=eye(3);
B2=[1 -1 0;0 0 -1;0 0 1];
B1=eye(3);

n= size(A,1);
m=size(B2,2);
  
S=[1 1 0;1 1 1;0 1 1]; %information structure

T=[1 1 0;1 1 1;0 0 1]; %chosen %
Rstruct=[1 1 0;1 1 0;0 0 1]; % Chosen R *** notice, cannot be computed with our simple algorithm
%% Performance
Q  =eye(n); R = eye(m);




[K0,J0_true,X0,Y0,J0_restriction] = StrucH2LMI_new_Gamma(A,B1,B2,Q,R,T,Rstruct,eye(n));      %Sparsity Invariance, obtain initial K0 to include later, Gamma=I
 %Lyapunov function
  
[Gamma,X,Y]  = compute_gamma(A,B1,B2,K0,T,Rstruct,Q,R);  %Compute Gamma

[K1,J1_true,X1,Y1,J1_restriction] = StrucH2LMI_new_Gamma(A,B1,B2,Q,R,T,Rstruct,Gamma);   % Now we use the Gamma we computed above.







