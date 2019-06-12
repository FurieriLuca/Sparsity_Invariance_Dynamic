function [ poly] = build_poly( cc )
s=sym('s');
n=length(cc);
poly=0*s;

for(i=1:n)
    poly=poly+s^(i-1)*cc(i);
end


end

