clear all clc

% Add Petsc "matlab-bin" directory path
addpath ~/Localdisk/petsc21april/matlab-bin/Matlab_Operators/PetscBinMatlab/


Title = 'Marmousi problem with optional damping';


% %% Choose problem
% freqV = [1 10 20 40];
% problem = menu('Choose a frequency: ', 'freq = 1', 'freq = 10', 'freq = 20', 'freq = 40');
% if ( problem > 0 )
%    freq = freqV(problem);
% end
% f = freq; 

%% Setting shifts in CSLP and Damping
f = 1; 
b1 = 1.0; b2 =1.0;                                  % Real & Imag. shifts in CSLP
alpha = 0.05;                                       % damping in Helmholtz eq.
gpwl = 10; 

prompt = {'frequency(1,10,20,40)','choose meshsize(10/20 gridpoint per wavelength)',...
'real shift in CSLP Preconditioner b1 = ','Imag. shift in CSLP preconditioner b2 = ',...
    'Damping parameter alpha = '};
defaults = {num2str(f), num2str(gpwl),num2str(b1), num2str(b2), num2str(alpha)};
numlines=1;
params=inputdlg(prompt,Title,numlines,defaults);

%parsing

f          = str2num(char(params(1)));
gpwl       = str2num(char(params(2)));
b1         = str2num(char(params(3)));
b2         = str2num(char(params(4)));
alpha      = str2num(char(params(5)));

%% Defining mesh sizes according to frequency

hV = [16 8 4 2]; 
if (f<=10)
h = hV(1);
elseif (f<=20)
h = hV(2);
else
h = hV(3); 
end

if (gpwl==20)
h = h/2;
end

 

%% Constructing discretized matrices for Helmholtz and CSLP and rhs vector. 

  [A,M,g,Nx,Ny] = buildMatrix(f,h,b1,b2,alpha);

 

% dx = 3000 - 952;                                    % size of domain x-axis 
% dy = 9200 - 1008;                                   % size od domain y-xais
% 
% hx = 16; hy = hx;                                   % step size in x- and y-direction

% nx = dx/hx; ny = dy/hy;                             % no of grid points in each direction

% nx = 2^9; ny = 2^11;
% hx = dx/nx; hy = dy/ny; 
 

%% adapt right-hand side vector. 
 
% Petsc would be configured with complex number, 
% it would not recognize real valued vectors or matrices,
% therefore make sure every matrix and vector passing to petsc
% must be complex valued. 

g(1) = g(1)+ 1e-12i; S{1}.b = g; 


%% Construct deflation vectors/Inter-grid transfer operators.

S{1}.A = A; 
S{1}.M = M; 

mynx = Nx; 
myny = Ny; 

nxvec = Nx; 
nyvec = Ny; 

l = 1; 

    S{l}.R = restrict2D(nxvec(l),nyvec(l));
    SizeVec_R(l) = length(S{l}.R);
    S{l}.R(1,1) =  S{l}.R(1,1)+ 1e-12i;
    
    S{l}.P = prolong2D(nxvec(l),nyvec(l));
    SizeVec_P(l) = length(S{l}.P);
    S{l}.P(1,1) =  S{l}.P(1,1)+ 1e-12i;

while mynx > 2 
   mynx = mynx/2; 
   myny = myny/2; 
   
   nxvec = [nxvec; mynx];
   nyvec = [nyvec; myny];
   
     l = l+1;
     
    S{l}.R = restrict2D(nxvec(l),nyvec(l));
    SizeVec_R(l) = length(S{l}.R);
    S{l}.R(1,1) =  S{l}.R(1,1)+ 1e-12i;
    
    S{l}.P = prolong2D(nxvec(l),nyvec(l));
    SizeVec_P(l) = length(S{l}.P);
    S{l}.P(1,1) =  S{l}.P(1,1)+ 1e-12i;
end

%% Saving to PETCs readable .DAT file 

% This will call PetscBinaryWrite.m function. Be sure to add path
% of the matlab-binary directory of Petsc

SizeVec_P(1) = SizeVec_P(1) + 1e-12i; 
NumLevel = SizeVec_P;

filename = ['f' num2str(f) 'gpWL' num2str(gpwl) 'a' num2str(alpha) '.dat'];

PetscBinaryWrite(filename,NumLevel, S{1}.A, S{1}.M, S{1}.P, S{1}.R, ...
            S{2}.P, S{2}.R, S{3}.P, S{3}.R, S{4}.P, S{4}.R, S{5}.P, ...
            S{5}.R, S{6}.P, S{6}.R, S{7}.P, S{7}.R, S{1}.b);
        fprintf('\n performing case when f = %d and alpha = %d \n',f,alpha);

 