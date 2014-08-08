
function g = rhs(N,Nx,Ny,hx,hy); 
%rhs_marm.m
% clear g; 
% Nx = nx; Ny = ny; 

 g = zeros(N,1);  
 g(round(N/2) - Nx/2)  = 1/(hx*hy); 
% g(round(N/2) + Nx/2 )  = 1/(hx*hy); 

% g = zeros(N,1);  
% g(end-round(Ny/2)) = 1/(hx*hy);




%  g(round(Nx/2)) = 1/(hx*hy);  

% g = zeros(N,1);  
% g(((Ny/2)