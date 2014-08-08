%  function grid = cnstructDefVec(nx,ny) 

nx = 752; 
ny = 2300; 

P1 = prolong2D(nx,ny);
R1 = restrict2D(nx,ny); 

% break; 

[nx2,ny2]= size(R1); 

P2 = prolong2D(nx2,ny2);
R2 = restrict2D(nx2,ny2);