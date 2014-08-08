function Z = restrict1D(n)

numelem = n;   % Number of elements 
Nh = numelem -1; NH = numelem/2-1;  % Number of fine and coarse grid points 
odd = [1:2:Nh]; even = [2:2:Nh-1]; % Indices of fine and coarse grid points 
                                   % on the fine grid  
Ne = length(odd);                 % Auliary data 

% Step 2: define deflation vectors W on 1D mesh 

u  = ones(Ne,1)/2;                 % Interpolation weights  
W  = zeros(Nh,NH);                 % Initialize deflation vectors  
W1 = eye(NH);                      % Interpolation to fine grid nodes lying 
                                   % on the coarse grid 
W2 = full(spdiags([u u],[-1 0],Ne,NH)); % Interpolation to fine grid nodes  
                                   % *not* lying on the coarse grid 
W(odd,:)  = W2;                    % Store 
W(even,:) = W1;                    % Idem 
W = (1/2)*sparse(W)'; 
vec1 = zeros(1,Nh);
W1 = [vec1;W];
vec2 = zeros(NH+1,1);
W2 = [vec2 W1];
vec3 = zeros(NH+1,1);
W3 = [W2 vec3];
vec4 = zeros(1,Nh+2);
W4 = [W3;vec4];
W = W4;
W(1,1) = 1;
W(NH+2,Nh+2) = 1;
W = W; 
Z = sparse(W);
