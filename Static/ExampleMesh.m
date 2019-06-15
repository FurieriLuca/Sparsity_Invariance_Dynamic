
%% Example 1: mesh network
%clc;
clear;close all

%Nn = [2,4,6,8,10];

Nn = [4];
L = 16; %number of agents with full information
%Nn = 4;

Num = length(Nn);
J   = zeros(Num,7);
degrees =zeros(1,2);
Index = Num;

%some inizialitizations
Jo1=0;
Jo2=0;
Jo1new=0;
Jo2new=0;
Jo1new2=0;
Jo2new2=0;
Jc=0;
leaders_vec=[14, 15, 3, 11, 2, 5, 9, 16, 8, 13, 7, 1, 12, 6, 10, 4];


n         = Nn(Index);
[Gp,Dist] = MeshGraph(n);  %% Plant Graph

%% Dynamics
A  = cell(n^2,n^2);   % matrices for A
B1 = cell(n^2,1);    % matrices for B

% Dynamical part
for i = 1:n^2
    B1{i} = [0;1];
    A{i,i} = [1 1; 1 2];
    for j = i+1:n^2
        if Gp(i,j) == 1
            A{i,j} = 15*exp(-norm(Dist(i,:)-Dist(j,:)))*eye(2); %%INCREASED THE COUPLINGS
            A{j,i} = 15*exp(-norm(Dist(i,:)-Dist(j,:)))*eye(2);
        end
    end
end
B2 = B1;

[As, B1s, B2s] = NetStateModel(A,B1,B2,Gp);


%%Initial communication graphs with 0 leaders
Gc = Gp;                             % mesh topology 
Gc2=MeshGraphReduced(n);             % maximal cliques



for leaders = 0 : 1 : L
    
    
    SP = bin(kron(Gc+eye(n^2),ones(1,2)));  %% sparsity patten +eye(n^2)
    SP2= bin(kron(Gc2+eye(n^2),ones(1,2)));
    
    %% Performance
    Q  = eye(2*n^2); R = eye(1*n^2);
    
    
    %% Structured Stablization P1
    %Ks = StruStaP1(A,B1,B2,Gp,Gc);    % controller
    %Js = sqrt(trace(lyap((As - B2s*Ks)',Q + Ks'*R*Ks)*(B1s*B1s'))); % H2 performance
    
    
    %% Structured Optimal control P2: LMI
    %Block-Diagonal Assumption
    [Ko1,Jo1,Jdiag] = StrucH2LMI(A,B1,B2,Gp,Q,R,SP);                     %Block Diagonal
    
    %Sparsity Invariance Approach
    [Ko1new,Jo1new,Jdiagnew,rnew] = StrucH2LMI_new(A,B1,B2,Gp,Q,R,SP);      %T=S  %rnew is the degree of separability of the Lyapunov function
    [Ko1new2,Jo1new2,Jdiagnew2,rnew2] = StrucH2LMI_new(A,B1,B2,Gp,Q,R,SP2); %T=T_new (maximal cliques) %rnew2 is the degree of separability of the Lyapunov function
    
    %% Structured Optimal control P2: gradient projection
    %[Ko2,Jo2,Iter] = StrucH2_Gradient(A,B1,B2,Gp,Q,R,Ko1,SP);                  %% Initial K from block-diagonal method
    %[Ko2new,Jo2new,Iter] = StrucH2_Gradient(A,B1,B2,Gp,Q,R,Ko1new,SP);  %% Initial K from sparsity invariance method
    
    %[Ko2new2,Jo2new2,Iter2] = StrucH2_Gradient(A,B1,B2,Gp,Q,R,Ko1new2,SP2);  %% Initial K from reduced cliques
    
    %% Augmented Lagragian Method by Lin et al
    %[Jaugl,Kaugl] = SH2_AugLag(As,B1s,B2s,Q,R,SP);
    %Jo3 = sqrt(Jaugl);
    
    %% LQR
    Kc = lqr(As,B2s,Q,R);
    Jc = sqrt(trace(lyap((As - B2s*Kc)',Q + Kc'*R*Kc)*(B1s*B1s')));
    
    
    %% decentralized control
    %Kd = DeceContr(A,B1,B2,Gp,1);
    %Jd = sqrt(trace(lyap((As - B2s*Kd)',Q + Kd'*R*Kd)*(B1s*B1s'))); % H2 performance
    
    
    %%randomly add a leader with full information to Gc and Gc2
    if(leaders<L)
       % [leaders_vec,r]=add_leader_random(leaders_vec,n^2)
        Gc(leaders_vec(leaders+1),:) = ones(1,n^2);
        Gc2(leaders_vec(leaders+1),:) = ones(1,n^2);
    end
    
    
    J(leaders+1,:) = [Jo1,Jo2,Jo1new,Jo2new,Jo1new2,Jo2new2,Jc];
    degrees(leaders+1,:)=[rnew,rnew2];
    
end

plots;






