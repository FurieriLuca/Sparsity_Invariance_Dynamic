equationZ=eye(m)+Y*Gz; %will force this to be stable

[NUM,DEN]=numden(equationZ);
maxdegnum=0;
for i=1:m
    for j=1:m
        if(size(coeffs(NUM(i,j),z),2)>maxdegnum)
            maxdegnum=size(coeffs(NUM(i,j),z),2);
        end
    end
end

freeZs = sym('freeZ',[m m*maxdegnum]);
freeZv = sdpvar(m,m*maxdegnum);   
for(i=1:maxdegnum)
    powers_z(i,1)=z^(i-1);
end



for i=1:m
    for j=1:m
        freepoly(i,j)=freeZs(i,[(j-1)*maxdegnum+1:j*maxdegnum])*powers_z;
    end
end


for i = 1:m
    for j = 1:m      
        fprintf('   Percentage %6.4f \n', 100*(m*(i-1)+j)/m/m );
        poles=roots(coeffs(DEN(i,j),z,'All'));
        U=1;
        for(k=1:size(poles))
            if(norm(poles(k))>1)
                U=U*(z-poles(k));
            end
        end
        polyno=NUM(i,j)-U*freepoly(i,j);
        cc=coeffs(polyno,z);
        [A_eq,b_eq]    = equationsToMatrix(cc,[vec(CQYs);vec(DQYs);vec(freeZs)]);      % Express system of equations in matrix form in terms of the vectorized versions of CQs and DQs
        A_eqs   = double(A_eq);    %A_eqs is the same as A_eq, for computation with sdpvars
        b_eqs    = double(b_eq);
        Constraints = [Constraints, A_eqs*[vec(CQYv);vec(DQYv);vec(freeZv)]== b_eqs]; 

    end
end
fprintf('Done \n')
clear freepoly
clear powers_s