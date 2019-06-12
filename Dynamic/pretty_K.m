[num,den]=numden(K);

%Compute the maximal order of any polynomial in K
maximal_coeff=0;
for(i=1:m)
    for(j=1:n)
        ccN_abs=double(abs(coeffs(num(i,j))));
        ccD_abs=double(abs(coeffs(den(i,j))));
        if(max(ccN_abs)>maximal_coeff)
            maximal_coeff=max(ccN_abs);
        end
        if(max(ccD_abs)>maximal_coeff)
            maximal_coeff=max(ccD_abs);
        end
    end
end

for(i=1:m)
    for(j=1:n)
        %ccN_rounded=round(double(coeffs(num(i,j)))/maximal_coeff,13);
        %ccD_rounded=round(double(coeffs(den(i,j)))/maximal_coeff,13);
        ccN_rounded=double(coeffs(num(i,j)))/maximal_coeff*100; %%rescales all coefficient to bring them into a reasonable range (same scaling for num and den, so the ratio does not change)
        ccD_rounded=double(coeffs(den(i,j)))/maximal_coeff*100;
        polyN=build_poly(ccN_rounded);
        polyD=build_poly(ccD_rounded);
        K_pretty(i,j)=simplifyFraction(polyN/polyD); %%Rewrites K into a nicer form
    end
end