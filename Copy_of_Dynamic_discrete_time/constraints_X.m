equationX=eye(n)+Gs*Y; %will force this to be stable
[NUM,DEN]=numden(equationX);
maxdegnum=0;
for i=1:n
    for j=1:n
        if(size(coeffs(NUM(i,j),s),2)>maxdegnum)
            maxdegnum=size(coeffs(NUM(i,j),s),2);  %This is the maximam degree of the numerator of the elements of I+GY
        end
    end
end

freeXs = sym('freeX',[n n*maxdegnum]); %%% These are the coefficients beta_{ijh}
freeXv = sdpvar(n,n*maxdegnum);   
for(i=1:maxdegnum)
    powers_s(i,1)=s^(i-1);
end



for i=1:n
    for j=1:n
        freepoly(i,j)=freeXs(i,[(j-1)*maxdegnum+1:j*maxdegnum])*powers_s;  %%defining free polynomials to factor out the unstable poles of I+GY 
    end
end


for i = 1:n
    for j = 1:n        
        fprintf('   Percentage %6.4f \n', 100*(n*(i-1)+j)/n/n );
        poles=roots(flip(coeffs(DEN(i,j),s)));
        U=1;
        for(k=1:size(poles))
            if(real(poles(k))>0)
                U=U*(s-poles(k));  %U contains the unstable poles
            end
        end
        polyno=NUM(i,j)-U*freepoly(i,j); %we require that NUM(i,j) has U(s) has a factor
        cc=coeffs(polyno,s);
        [A_eq,b_eq]    = equationsToMatrix(cc,[vec(CQYs);vec(DQYs);vec(freeXs)]);      
        A_eqs   = double(A_eq);    
        b_eqs    = double(b_eq);
        Constraints = [Constraints, A_eqs*[vec(CQYv);vec(DQYv);vec(freeXv)]== b_eqs];        
    end
end
fprintf('Done \n')
clear freepoly
clear powers_s