% vel_marm.m   for h = 4 

function get_output = vel_marm4h(Nx,Ny)
 tic
% clear all 
addpath deflation_vectors/

fid  = fopen('velocity.h@','r','l');
vel_vec = fread(fid,'single');
vel_mat = reshape(vel_vec,751,2301);
% trim_vel_mat = vel_mat(239:751, 253:2301); 
trim_vel_mat = vel_mat(1:751-238, 1:2301-252); 
% Trimed from west and north to make it in-powers-of 2. 
    % Nx = 2^9 and Ny = 2^11
trim_vel_vec = reshape(trim_vel_mat,513*2049,1); 
fclose(fid); 

NxOrig = 2^9; NyOrig = 2^11; 
% 
% R = restrict2D(NxOrig,NyOrig);
% Nxh8 = NxOrig/2; Nyh8 = NyOrig/2; 
% trim_vel_vec_h8 = R*trim_vel_vec; clear R;
% 
% trim_vel_mat_h8 = reshape(trim_vel_vec_h8,Nxh8 +1,Nyh8 +1); 

trim_vel_mat_required = trim_vel_mat; 
trim_vel_vec_required = trim_vel_vec; 


% figure, imagesc(trim_vel_mat_required), colorbar
kold = trim_vel_mat_required;
dif = max(max(trim_vel_mat)) - min(min(trim_vel_mat));
knew = dif + 0.5 * (kold - dif); 
% figure, imagesc(knew), colorbar
knew1 = dif + 0.25 * (kold - dif); 
% figure, imagesc(knew1), colorbar


get_output = reshape(knew1,(Nx +1)*(Ny +1),1); 

