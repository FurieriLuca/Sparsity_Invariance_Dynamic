
% This is a function copied from Lin et al's wrok

% Created 2010-May-05
% alternating iterative method between to minimize the augmented Lagrangian
% inputs: initial step F, dual variable V, system data, A,B1,B2,C2,Q,R
% penalty weight c, complementary structural identity and the tolerance
% outputs: minimizer F, optimal augmented Lagrangian value Lag, iteration
% numbers iterNo.

function [F,Lag,iterNo] = SH2_AugLag_ALT(Fold,V,A,B1,B2,C2,Q,R,c,ISc,tol)

% backtracking line search data

    alpha = 0.3; 
    beta = 0.5;

% data setup
    
    iterNo = 0;
    LagEvalNo=0;
    maxIter = 100;

    F = Fold;

%check if F is a stabilizing feedback gain

    maxEigAcl  = max(real(eig(A-B2*F*C2)));
    if maxEigAcl > 0
        error('input F is not stabilizing!')
    end


while 1 
    
    L = lyap(A-B2*F*C2,B1*B1');   
    P = lyap((A-B2*F*C2)',Q + C2'*F'*R*F*C2);
    Lag = trace( L*(Q+C2'*F'*R*F*C2) ) + trace(F'*V)+ 0.5*c*norm(F.*ISc,'fro')^2;
    LagEvalNo = LagEvalNo+1;
    gradAugLag = 2*(R*F*C2-B2'*P)*L*C2' + V + c*(F.*ISc);    
    
    % if the number of alternating iterative steps is larger than the
    % maximum number, the algorithm terminates and return the current value
    % of augmented Lagrangian
    
    ngrad = norm(gradAugLag,'fro');
    
    if iterNo>=maxIter
        disp(['maximum number ' num2str(maxIter) ' of Augmented Lagrangian iterations reached! The norm of gradient is ',...
            num2str(ngrad,'%11.4E'),' Lagrangian is ',num2str(Lag,'%11.4E')]);
        break;
    end
    
% %     ------------------------------------------------------------------
% %     ---------------Solve vectorized linear system-------------
    
% %     % compute the solution of the F equation depending on L, P and then
% %     % compute the difference of this new F and the previous F to get the 
% %     % direction Ftilda

    Ftilda = SH2_AugLag_ALT_solveF(P,L,V,B2,C2,R,ISc,c) - F;      

% %     ------------------------------------------------------------------
% %     ---------------Conjugate gradient method-------------

%     Ft0 = zeros(size(F));
%     Ftilda = SH2_AugLag_ALT_CG(R,C2,L,gradAugLag,ISc,Ft0,c);
    
    % compare the solution of CG method with direction method
%     norm(Ftilda1-Ftilda)

% %     ------------------------------------------------------------------

    iterNo = iterNo+1;
    
    % check if Ftilda is a descent direction;
    
    if trace(Ftilda'*gradAugLag) >= 0 
        error('dirF is not a descent direction!')
    end
    
    % stopping criterion is the inner product between Ftilda and the
    % gradient is small enough.

    nFtildaGrad = abs(trace(Ftilda'*gradAugLag));
    
    if nFtildaGrad<tol
        break;
    end  
    
    % use backtracking line search and also check closed-loop stability
    
    stepsize = 1;
    
    while 1
        
        Ftemp = F + stepsize*Ftilda;
        maxEigAcltemp  = max(real(eig(A-B2*Ftemp*C2)));
        Ltemp = lyap(A-B2*Ftemp*C2,B1*B1');
        Lagtemp = trace( Ltemp*(Q+C2'*Ftemp'*R*Ftemp*C2) ) + trace(Ftemp'*V)+ 0.5*c*norm(Ftemp.*ISc,'fro')^2;
        LagEvalNo = LagEvalNo+1;
        if maxEigAcltemp < 0 && Lag - Lagtemp > stepsize*alpha*trace(-Ftilda'*gradAugLag)
            break;
        end
        
        stepsize = stepsize*beta;
        
        if stepsize<1.e-16
            maxEigAcltemp %#ok<NOPRT>
            error('extremly small stepsize in alternating iterative method in minimizing the Augmented Lagrangian!');            
        end        
    end 
    
    % update current step
    
    F = Ftemp;
    Lag = Lagtemp;
end
