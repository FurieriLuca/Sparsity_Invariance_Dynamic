% Created 2010-May-05
% Newton-CG method to minimize the augmented Lagrangian, 
% L = J(K) + trace(V'*K) + 0.5*c*|| K.*ISc ||

function [K,g,iterNo,LagEvalNo,totalCGstep] = SH2_AugLag_newton(Kold,V,A,B1,B2,C2,Q,R,ISc,c,tol)

% backtracking data

    alpha = 0.3;
    beta = 0.3;   

% evaluate augmented Lagrangian    
    
    K = Kold;
    L = lyap(A-B2*K*C2,B1*B1');
    Lag = trace( L*(Q+C2'*K'*R*K*C2) ) + trace(K'*V) + 0.5*c*norm(K.*ISc,'fro')^2;
    LagEvalNo = 1;
    
% maximum number of iterations of newton steps

    maxIter = 10;  
    totalCGstep=0;
    iterNo = 0;

while 1
    
    % compute the gradient direction
    
    L = lyap(A-B2*K*C2,B1*B1');               
    P = lyap((A-B2*K*C2)',Q + C2'*K'*R*K*C2); 
    Kgrad = V + 2*(R*K*C2-B2'*P)*L*C2' + c*K.*ISc;
    nKgrad = norm(Kgrad,'fro');
    
    % if gradient is small enough, stop
    
    if nKgrad < tol 
        break;
    end    
    
    % compute the Newton direction
    
    % always use zero initial condition for CG !!
    
    Kt0 = zeros(size(K));   
    [Knt,res,CGstep] = SH2_AugLag_newton_CG(A,B2,C2,K,P,L,R,Kgrad,ISc,Kt0,c);     

        
    totalCGstep = totalCGstep+CGstep;
    
    % if the number of Newton steps is larger the maximum
    % terminate the algorithm and give the best optimal value
        
    if iterNo>=maxIter
        disp(['maximum number ' num2str(maxIter) ' of newton iterations reached! The norm of gradient is ',...
            num2str(nKgrad,'%11.4E'),' Lagrangian is ',num2str(Lagtemp,'%11.4E')]);
        break;
    end
    
    stepsize = 1;
    iterNo = iterNo + 1;
    
    % use backtracking line search 
    
    while 1
        
        Ktemp = K + stepsize*Knt;
        maxEigAcltemp  = max(real(eig(A-B2*Ktemp*C2)));
        Ltemp = lyap(A-B2*Ktemp*C2,B1*B1');
        Lagtemp = trace( Ltemp*(Q+C2'*Ktemp'*R*Ktemp*C2) ) + trace(Ktemp'*V)+ 0.5*c*norm(Ktemp.*ISc,'fro')^2;
        LagEvalNo = LagEvalNo+1;
        if maxEigAcltemp<0 && Lag - Lagtemp > stepsize*alpha*trace(-Knt'*Kgrad)
            break;
        end
        stepsize = stepsize*beta;
        if stepsize<1.e-16
            maxEigAcltemp %#ok<NOPRT>
            error('extremly small stepsize in newton method for Lagrangian!');            
        end
    end
    
    % update current step
    
    K = Ktemp; Lag = Lagtemp;
end

% dual objective function

g = trace(L*(Q + C2'*K'*R*K*C2)) + trace(K'*V) + 0.5*c*norm(K.*ISc,'fro')^2;  
