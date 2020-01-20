equationW=(eye(n)+Gs*Y)*Gs; %will force this to be stable

[NUM,DEN]=numden(equationW);
maxdegnum=0;
for i=1:n
    for j=1:m
        if(size(coeffs(NUM(i,j),s),2)>maxdegnum)
            maxdegnum=size(coeffs(NUM(i,j),s),2);
        end
    end
end

freeWs = sym('freeW',[n m*maxdegnum]);
freeWv = sdpvar(n,m*maxdegnum);   
for(i=1:maxdegnum)
    powers_s(i,1)=s^(i-1);
end



for i=1:n
    for j=1:m
        freepoly(i,j)=freeWs(i,[(j-1)*maxdegnum+1:j*maxdegnum])*powers_s;
    end
end


for i = 1:n
    for j = 1:m      
        fprintf('   Percentage %6.4f \n', 100*(m*(i-1)+j)/m/n );
        poles=roots(flip(coeffs(DEN(i,j),s)));
        U=1;
        for(k=1:size(poles))
            if(real(poles(k))>0)
                U=U*(s-poles(k));
            end
        end
        polyno=NUM(i,j)-U*freepoly(i,j);
        cc=coeffs(polyno,s);
        [A_eq,b_eq]    = equationsToMatrix(cc,[vec(CQYs);vec(DQYs);vec(freeWs)]);      % Express system of equations in matrix form in terms of the vectorized versions of CQs and DQs
        A_eqs   = double(A_eq);    %A_eqs is the same as A_eq, for computation with sdpvars
        b_eqs    = double(b_eq);
        Constraints = [Constraints, A_eqs*[vec(CQYv);vec(DQYv);vec(freeWv)]== b_eqs]; 

    end
end
fprintf('Done \n')
clear freepoly
clear powers_s