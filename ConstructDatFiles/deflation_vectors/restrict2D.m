% Builds the 2D RESTRICTION operator,
% Code can be converted for 2D using a kronecker delta product of 
% the Restriction matrix on a 1D regular grid
% IT IS SUPPOSED THAT BOUNDARY NODES ARE INCLUDED
% IN DISCRETIZATION SCHEME. HENCE THE TRANSFORMATION WILL BE IDENTICAL TO
% IDENTITY ON BOUNDARY NODES.
% bilinear interpolation is used 
% STANDARD WEIGHTAGE IS USED.
%                                                           3 MARCH 2011
% ======================================================================
function Z =  restrict2D(nx,ny)          % n is stepsize s.t. h=1/n


Zx = restrict1D(nx);
Zy = restrict1D(ny);
Z = kron(Zy,Zx); 
ZZ = sparse(Z);

