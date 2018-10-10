function [Gc] = MeshGraphReduced(n,L)
%This function creates a graph with  only the maximum cliques (dim 2) of the mesh
%graph
if(mod(n,2)==1)
        disp('here we suppose n is an even number for simplicity. Can be extended')
end

%Dist(:,1) = X(:)*1.2; Dist(:,2) = Y(:)*1.2;  %% coordinates

Gc=kron(eye(n^2/2),ones(2,2));
for(i=1:L)
        Gc(i,:) = ones(1,n*n);
end

end

