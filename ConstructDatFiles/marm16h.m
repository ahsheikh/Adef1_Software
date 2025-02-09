%*****************  -DELTA(u) -(b1 - i*b2)*k(x) u = g  ********** 
%***************    f = deltafunction(x-300,y) ***********
%************** DOMAIN: (0,600)X(0,1000)  ********
%************ DOMAIN is bounded by first order radiation Bd.Cond.
% this code gives WEDGE helm. op on a rectangular grid and 
% works fine, as the solution so computed from discretized
% operator and exact solution are same.   14 dec 2010



function [helm,cslp,trim_vel_vec, N,hx,hy] = marm16h(f,hx,hy,b1,b2,alpha)
%  tic
% clear all,  close all, clc

% f = 40; 
% Nx = 2^9; Ny = 2^11; 

% Nx = 750; Ny = 2300;
dx = 3000 - 952; dy = 9200-1008;
% Nx = 115; Ny = 190;
Nx = dx/hx; 
Ny = dy/hy; 
% hx = dx/Nx; 
% hy = dy/Ny;

N = (Nx+1)*(Ny+1); 

A = sparse(N,N);
% Id = speye(N);
% ------ CALLING velocity  OVER GRID -----
% d1=400; d2=500; d3=800; d4=600;
% vel_marm; 
[trim_vel_vec] = vel_marm16h(Nx,Ny); 
kmat = trim_vel_vec;  
kmatt = (2*f*pi)./kmat; %kmatt = (b1-1i*b2)*kmatt;
% -----------------------------------------------------
% sizekmatt = length(kmatt);
% A1 = ((2*(hx^2)) + (2*(hy^2)))*(speye(N)) - ...
% ( (b1-1i*b2) )*( (hx^2)*(hy^2) )*( sparse(1:sizekmatt,1:sizekmatt,kmatt.^2,sizekmatt,sizekmatt));%  diag(kmatt.^2)));
% A = A + A1; clear A1;

E = sparse(1:N-Nx-1,Nx+2:N,- hx^2,N,N);
EE = sparse(Nx+2:N,1:N-Nx-1,- hx^2,N,N);
A = A+ E+EE;

clear E, clear EE;                                                                                                                                           

E = sparse(2:N,1:N-1,- hy^2,N,N);                                                                                                                            
EE= sparse(1:N-1, 2:N, - hy^2, N,N);                                                                                                                         
A =A+ E + EE;                                                                                                                                                
clear E, clear EE;                                                                                                                                           

%              lower boundary                                                                                                                                
E = sparse(1:Nx+1, 1:Nx+1,2*hy*(hx^2)*sqrt(-1)*kmatt(1:Nx+1), N,N);                                                                                          
EE = sparse(1:Nx+1, Nx+1+1:Nx+1+Nx+1, - hx^2, N,N);                                                                                                          
A = A+ E + EE;                                                                                                                                               
clear E, clear EE;                                                                                                                                           

%              right boundary                                                                                                                               
E = sparse(Nx+1:Nx+1:N, Nx+1:Nx+1:N,2*hx*(hy^2)*sqrt(-1)*kmatt(Nx+1:Nx+1:N),N,N);                                                                           
EE = sparse(Nx+1:Nx+1:N,Nx:Nx+1:N-1, - hy^2, N,N);                                                                                                          
A = A + E + EE;                                                                                                                                             
clear E; clear EE;                                                                                                                                          

%              upper boundry                                                                                                                                
E = sparse( (Nx+1)*Ny+1:(Nx+1)*Ny+Nx+1, (Nx+1)*Ny+1:(Nx+1)*Ny+Nx+1,2*hy*sqrt(-1)*(hx^2)*kmatt( (Nx+1)*Ny+1:(Nx+1)*Ny+Nx+1),N,N);                            
EE= sparse( (Nx+1)*Ny+1:(Nx+1)*Ny+(Nx+1), (Nx+1)*(Ny-1) + 1: (Nx+1)*(Ny-1)+Nx+1, - hx^2, N,N);                                                              
A = A + E + EE;                                                                                                                                             
clear E; clear EE;                                                                                                                                          
%             left boundary                                                                                                                                
E = sparse(1:Nx+1:N-1,1:Nx+1:N-1, 2*hx*(hy^2)*sqrt(-1)*kmatt(1:Nx+1:N-1),N,N);                                                                             
EE= sparse(1:Nx+1:N-1,2:Nx+1:N, - hy^2, N,N);                                                                                                              
A = A + E + EE;                                                                                                                                            
clear E, clear EE;

%% Inserting speed-profile in Helmholtz with damping. 

sizekmatt = length(kmatt);
A1_helm = ((2*(hx^2)) + (2*(hy^2)))*(speye(N)) - ...
( (b1-1i*b2) )*( (hx^2)*(hy^2) )*( sparse(1:sizekmatt,1:sizekmatt,kmatt.^2,sizekmatt,sizekmatt));%  diag(kmatt.^2)));
helm = A + A1_helm; clear A1;

%% Inserting speed-profile in CSLP Preconditioner. No Damping.

A1_cslp = ((2*(hx^2)) + (2*(hy^2)))*(speye(N)) - ...
( (1 - alpha) )*( (hx^2)*(hy^2) )*( sparse(1:sizekmatt,1:sizekmatt,kmatt.^2,sizekmatt,sizekmatt));%  diag(kmatt.^2)));

cslp = A + A1_cslp ; 
%% Vanishing points alongwith boundaries. 

%A(Nx+1:Nx+1:N-1, Nx+2:Nx+1:N) = 0;                                                                                                                         
%A(Nx+2:Nx+1:N, Nx+1:Nx+1:N-1) = 0;  

% toc, tic
% matlabpool(4);              % starts parallel in 2 processors
% parfor j = Nx+1:Nx+1:N-1
%     A(j,j+1) = 0;
%     A(j+1,j) = 0;
% end
% matlabpool close; 

% toc, 
tic

for j = Nx+1:Nx+1:N-1
    
    helm(j,j+1) = 0;                % for Helmholtz
    helm(j+1,j) = 0;
    
    cslp(j,j+1) = 0;                % for CSLP
    cslp(j+1,j) = 0;
end

toc

% for j = Nx+1:Nx+1:N-1
%     A(j,j+1) = 0;
%     A(j+1,j) = 0;
% end

helm=(1/((hx^2)*(hy^2)))*helm;
cslp=(1/((hx^2)*(hy^2)))*cslp;
% A=(1/((hx^2)*(hy^2)))*A;
