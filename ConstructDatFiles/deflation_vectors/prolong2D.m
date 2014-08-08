% Builds the 2D Prolongation operator with different nx and ny nodes
% Code can be converted for 2D using a kronecker delta product of 
% the Prolongation matrix on a 1D regular grid
% IT IS SUPPOSED THAT BOUNDARY NODES ARE INCLUDED
% IN DISCRETIZATION SCHEME. HENCE THE TRANSFORMATION WILL BE IDENTICAL TO
% IDENTITY ON BOUNDARY NODES.
% bilinear interpolation is used 
%                                                           3 MARCH 2011
% ======================================================================

function Z =  prolong2D(nx,ny)

Zx = prolong1D(nx);  [nxtmp,nytmp] = size(Zx); nytmp;
Zy = prolong1D(ny);  [nxtmp,nytmp] = size(Zy); nytmp; 
Z = kron(Zy,Zx); 
ZZ = sparse(Z);

% If desired, do some checks 1D. Anything goes here. 

if (0)

 figure
 imagesc(Zx), colorbar 
 xlabel('Number of coarse grid point')
 ylabel('Number of fine grid point')
 title('Multigrid deflation vectors on 1D grid')

 figure
 imagesc(Zx), colorbar 
 xlabel('Number of coarse grid point')
 ylabel('Number of fine grid point')
 title('Multigrid deflation vectors on 1D grid')

end




% % If desired, do some checks 2D. Anything goes here
if (0)
 figure
 imagesc(Z), colorbar
 xlabel('Number of coarse grid point')
 ylabel('Number of fine grid point')
 title('Multigrid deflation vectors on 2D grid')
end 






