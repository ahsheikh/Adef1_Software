function ZZ =  mg_def_vec(nx,ny)
Wx = mg_def(nx);
Wy = mg_def(ny);
Z = kron(sparse(Wy),sparse(Wx)); 
% Z = kron2(Wy,Wx); 
ZZ = sparse(Z);

