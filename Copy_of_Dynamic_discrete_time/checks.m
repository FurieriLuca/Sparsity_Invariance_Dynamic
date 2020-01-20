
%CHECKS


%{
disp('Performing the check of stability of I+Gs*Y \n')
for i = 1:n
    for j = 1:n
        fprintf('   Percentage %6.4f \n', 100*(n*(i-1)+j)/n/n );
        if(i==3&&j==5)
            disp('stop debugger')
        end
        [num,den]=numden(eq(i,j));
        poles=roots(flip(coeffs(den,s)));
        zeros=roots(flip(coeffs(num,s)));
        num_nogain=1;
        den_nogain=1;
        for(k=1:size(zeros))
            num_nogain=num_nogain*(s-zeros(k));
        end
        for(k=1:size(poles))
            den_nogain=den_nogain*(s-poles(k));
        end
        gain_num=simplify(num/num_nogain);
        gain_den=simplify(den/den_nogain);
        
        zeros_rounded=round(double(zeros),4);
        poles_rounded=round(double(poles),4);
        num_nogain_rounded=1;
        den_nogain_rounded=1;
        for(k=1:size(zeros))
            num_nogain_rounded=num_nogain_rounded*(s-zeros_rounded(k));
        end
        for(k=1:size(poles))
            den_nogain_rounded=den_nogain_rounded*(s-poles_rounded(k));
        end
        num_rounded=num_nogain_rounded*gain_num;
        den_rounded=den_nogain_rounded*gain_den;
        Xr(i,j)=num_rounded/den_rounded;
    end
end
%}


%Check stability parameters.  
tollerance=1e-3;

eq1 = inv(eye(n)-Gz*Yr*inv(eye(n)+Gz*Yr)); %% This is (I-G*K)^-1=(I-G*Y*X^-1)^-1=(I-G*Y*(I+G*Y)^-1)^-1
for i = 1:n
    for j = 1:n
        fprintf('   Percentage %6.4f \n', 100*(n*(i-1)+j)/n/n )
        
        if(isstable(syms2tf(eq1(i,j)))==0)   %% if it says it's unstable... check that there are very close zero/pole cancellations, so it is actually stable
            fprintf('   Apparently unstable: check close zero/pole cancellations\n')
            [num,den]=numden(eq1(i,j));
            poles=double(roots(flip(coeffs(den,s))));
            zeroes=double(roots(flip(coeffs(num,s))));
            
            unstable_poles=[];
            for(k=1:size(poles))
                if(real(poles(k))>=0)
                    unstable_poles=[unstable_poles, poles(k)];
                end
            end
            stable=zeros(1,size(unstable_poles,2));
            zeroes_comparison=zeroes;
            for(l=1:size(unstable_poles,2))
                done=0;
                for(v=1:size(zeroes,1))
                    if(abs(zeroes_comparison(v)-unstable_poles(l))<tollerance &&done==0) 
                        stable(l)=1;
                        zeroes_comparison(v)=inf; %not to use it again for the next unstable poles, this zero is already used.
                        done=1;
                    end
                end
            end
            for(z=1:size(stable,2))
                if(stable(z)==0)
                    disp('******ACTUALLY UNSTABLE :-( ******')
                end
            end
            
        end
    end
end


eq2 = Yr*inv(eye(n)+Gs*Yr)*inv(eye(n)-Gs*Yr*inv(eye(n)+Gs*Yr));  %This is the check for Y=K(I-GK)^{-1}%
for i = 1:m
    for j = 1:n
        fprintf('   Percentage %6.4f \n', 100*(n*(i-1)+j)/n/n )
        
        if(isstable(syms2tf(eq2(i,j)))==0) %% if it says it's unstable... check that there are very close zero/pole cancellations, so it is actually stable
            fprintf('   Apparently unstable: check close zero/pole cancellations\n')
            [num,den]=numden(eq2(i,j));
            poles=double(roots(flip(coeffs(den,s))));
            zeroes=double(roots(flip(coeffs(num,s))));
            
            unstable_poles=[];
            for(k=1:size(poles))
                if(real(poles(k))>=0)
                    unstable_poles=[unstable_poles, poles(k)];
                end
            end
            stable=zeros(1,size(unstable_poles,2));
            zeroes_comparison=zeroes;
            for(l=1:size(unstable_poles,2))
                done=0;
                for(v=1:size(zeroes,1))
                    if(abs(zeroes_comparison(v)-unstable_poles(l))<tollerance && done==0)
                        stable(l)=1;
                        zeroes_comparison(v)=inf; 
                        done=1;
                    end
                end
            end
            for(z=1:size(stable,2))
                if(stable(z)==0)
                    disp('******ACTUALLY UNSTABLE :-( ******')
                end
            end
            
        end
    end
end


eq3 = inv(eye(n)-Gz*Yr*inv(eye(n)+Gz*Yr))*Gz;   %%This is the check for (I-GK)^{-1}G
for i = 1:n
    for j = 1:m
        fprintf('   Percentage %6.4f \n', 100*(n*(i-1)+j)/n/n )
        
        if(i==2 && j==2)
            disp('ciao')
        end
        if(isstable(syms2tf(eq3(i,j)))==0) %% if it says it's unstable... check that there are very close zero/pole cancellations, so it is actually stable
            fprintf('   Apparently unstable: check close zero/pole cancellations\n')
            [num,den]=numden(eq3(i,j));
            poles=double(roots(flip(coeffs(den,s))));
            zeroes=double(roots(flip(coeffs(num,s))));
            
            unstable_poles=[];
            for(k=1:size(poles))
                if(real(poles(k))>=0)
                    unstable_poles=[unstable_poles, poles(k)];  
                end
            end
            stable=zeros(1,size(unstable_poles,2));
            zeroes_comparison=zeroes;
            for(l=1:size(unstable_poles,2))
                done=0;
                for(v=1:size(zeroes,1))
                    if(abs(zeroes_comparison(v)-unstable_poles(l))<tollerance && done==0) 
                        stable(l)=1;
                        zeroes_comparison(v)=inf; 
                        done=1;
                    end
                end
            end
            for(z=1:size(stable,2))
                if(stable(z)==0)
                    i
                    j
                    disp('******ACTUALLY UNSTABLE :-( ******')
                end
            end
            
        end
    end
end


eq4 = inv(eye(n)-Yr*inv(eye(n)+Gs*Yr)*Gs);  %This is the check for (I-KG)^{-1}
for i = 1:m
    for j = 1:m
        fprintf('   Percentage %6.4f \n', 100*(n*(i-1)+j)/n/n ) 
        if(isstable(syms2tf(eq4(i,j)))==0)   %% if it says it's unstable... check that there are very close zero/pole cancellations, so it is actually stable
            fprintf('   Apparently unstable: check close zero/pole cancellations\n')
            [num,den]=numden(eq4(i,j));
            poles=double(roots(flip(coeffs(den,s))));
            zeroes=double(roots(flip(coeffs(num,s))));
            
            unstable_poles=[];
            for(k=1:size(poles))
                if(real(poles(k))>=0)
                    unstable_poles=[unstable_poles, poles(k)];
                end
            end
            stable=zeros(1,size(unstable_poles,2));
            zeroes_comparison=zeroes;
            for(l=1:size(unstable_poles,2))
                done=0;
                for(v=1:size(zeroes,1))
                    if(abs(zeroes_comparison(v)-unstable_poles(l))<tollerance && done==0) 
                        stable(l)=1;
                        zeroes_comparison(v)=inf;
                        done=1;
                    end
                end
            end
            for(z=1:size(stable,2))
                if(stable(z)==0)
                    disp('******ACTUALLY UNSTABLE :-( ******')
                end
            end
            
        end
    end
end


