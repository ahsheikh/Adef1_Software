function Z=  prolong1D(no)
numelem = no;   % Number of elements 
Nh = numelem+1; NH = numelem/2+1;  % Number of fine and coarse grid points 
odd = [1:2:Nh]; even = [2:2:Nh-1]; % Indices of fine and coarse grid points 
                                   % on the fine grid  
Ne = length(even);                 % Auliary data 

% Step 2: define deflation vectors W on 1D mesh 

u  = ones(Ne,1)/2;                 % Interpolation weights  
W  = zeros(Nh,NH);                 % Initialize deflation vectors  
W1 = eye(NH);                      % Interpolation to fine grid nodes lying 
                                   % on the coarse grid 
W2 = full(spdiags([u u],[0 1],Ne,NH)); % Interpolation to fine grid nodes 
                                   % *not* lying on the coarse grid 
W(odd,:)  = W1;                    % Store 
W(even,:) = W2;                    % Idem 
W = sparse(W);
Z = W; 