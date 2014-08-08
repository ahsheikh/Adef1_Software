% buildMatrix.m
function [A,M,g,Nx,Ny] = buildMatrix(f,h,b1,b2,alpha);


dx = 3000 - 952;                                    % size of domain x-axis 
dy = 9200 - 1008;

hx = h;                     % equidistant grids in x- and y-directions
hy = h; 


%%  Generate matrices for Helmholtz and CSLP

switch h
case 16
[A,M,trim_vel_vec, N,hx,hy] = marm16h(f,hx,hy,b1,b2,alpha);
case 8 
[A,M,trim_vel_vec, N,hx,hy] = marm8h(f,hx,hy,b1,b2,alpha);
case 4
[A,M,trim_vel_vec, N,hx,hy] = marm4h(f,hx,hy,b1,b2,alpha);
case 2
[A,M,trim_vel_vec, N,hx,hy] = marm2h(f,hx,hy,b1,b2,alpha);
end

%% Construct rhs ; point source. 
Nx = dx/hx; Ny = dy/hy; 

g = rhs(N,Nx,Ny,hx,hy); 

  
          
          