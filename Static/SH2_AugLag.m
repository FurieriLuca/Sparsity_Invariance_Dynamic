function [Jopt,K] = SH2_AugLag(A,B1,B2,Q,R,IS)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% This function is copied from Fu Lin et al's work

% Created 2010-May-05
% the augmented Lagrangian method to solve the Structured H2 problem
% inputs: system data A,B1,B2,C2,Q,R,structural identity IS,
% outputs: optimal structured gain Kopt, optimal dual variable Vopt, and
% Jopt

% the code minimizes the augmented Lagrangian 
% L = J(K) + trace(Enu'*K) + 0.5*c*||K.*ISc||^2
% the last term is the quadratic penalty on the structured constraint that
% h(K) = K.*ISc = 0
% For every dual variable Enu, the augmented Lagrangian is minimized by
% either Newton's method or the alternating iterative method,
 

% record the number of CG steps, NT steps, and the alternating iterative
% steps

    grandCGstep = 0;
    ntIterNo = 0;
    
    ALTIterNo = 0;

% parameters for the update of the weight c
% c+ = gam*c
% and c stops increasing if c >= tau

    c = 1;
    gam = 5;
    tau = 1.e+5;

% complementary set of the constraint

    ISc = ones(size(IS)) - IS;   

% the tolerance of the stopping criterion

    tol = 1.e-5;
    
% parameters of backtracking line search

    beta = 0.1;
    alpha = 0.3;
    
    C2 = eye(size(A));    %% For our example, all the simulation is state feedback, by Yang

% check the number of output 

    [pn,nn] = size(C2);

% case 1, if C2 is skinny full column-rank, initialize the descent method with 
% K = Kc*inv(C2'*C2)*C2';

if pn >= nn
    
    % compute the centralized gain    
    
     [Kc,S,E] = lqr(A,B2,Q,R);  %% This is not efficient and not scalable.  by Yang
     K0 = Kc*inv(C2'*C2)*C2'; 
end

% case 2, if C2 is fat full row-rank, do the coordinate transformation

if pn < nn
    
    T = [C2; null(C2)'];
    invT = inv(T);
    A = T*A*invT;  
    B1=T*B1; B2=T*B2;
    Q = invT'*Q*invT;  
    [Kc,S,E] = lqr(A,B2,Q,R);

    % after the change of coordiantes, C2 = [I O];
    % then C2 should be eliminated, but to make the code
    % consistent we have C2=I
    % the initial condition is the optimal feedback gain Kc
    
    C2 = eye(size(A));
    IS = [IS, zeros(pn,nn-pn)];    
    ISc = ones(size(IS)) - IS;
    K0 = Kc;
end
       
    K = K0;
    stepno = 0;

% start dual variable V with zeros    
    
    V = zeros(size(K));  

% gradient direction of the dual function    

    dV = K.* ISc;  
    ndV = norm(dV,'fro');

    disp('stepno      penaltyWeight         normGradDualVal           dualFunVal        ')
    
% tic
while 1
    
    % the stepsize for updating V
    % V+ = V + stepsize*dV

    stepsize = c;    
    
    V_temp = V + stepsize*dV;       
    
    % g is the dual function, i.e., minimum of augmented Lagrangian
    
    % % ----------------------------------------------------------------
    % % ---------------------Newton's method---------------------------    
    
     [Ktemp,gtemp,iterNo,LagEvalNo,totalCGstep] = SH2_AugLag_newton(K,V_temp,A,B1,B2,C2,Q,R,ISc,c,tol);  
     grandCGstep = grandCGstep + totalCGstep;
     ntIterNo = ntIterNo + iterNo ;
    
    % % ----------------------------------------------------------------    
    % % ---------------------Alternating iterative method---------------
    
     %[Ktemp,gtemp,iterNo] = SH2_AugLag_ALT(K,V_temp,A,B1,B2,C2,Q,R,c,ISc,tol);     
     %ALTIterNo = ALTIterNo + iterNo;

    % % ----------------------------------------------------------------    

    % increase c by ratio gamma until it is large enough, i.e., c>=tau

    if c <= tau
        c=c*gam;
    end

%     totalLagEvalNo = totalLagEvalNo + LagEvalNo;
%     if gtemp>Jopt
%         disp('dual function is greater than local minimum')
%     end    
    
    % update the current step
    
    V = V_temp; 
    g = gtemp;
    K = Ktemp;
    stepno = stepno + 1;
    dV = K.* ISc;  % gradient direction of Enu
    ndV = norm(dV,'fro');
    
    % display the norm of dEnu, the dual function value and the local optimal
    
    disp([num2str(stepno),'           ',num2str(c,'%11.1E'),'              ',...
        num2str(ndV,'%11.4E'),'             ', num2str(g,'%11.4E')])
    
    % stopping criterion is the gradient of the dual variable small enough
    
    if ndV < tol
        break;
    end
end
% toc

Jopt = g;


end

