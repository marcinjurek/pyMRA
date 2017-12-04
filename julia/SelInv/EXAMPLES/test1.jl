nnodes = 8;
nnodes = convert(Int32,nnodes);
nnz = 12;
testM = eye(nnodes)*4;
testM[2,1] = 1
testM[1,2] = 1
testM[3,2] = 1.2
testM[2,3] = 1.2
testM_sparse = sparse(testM);

function SelInv(testM_sparse)

colptr = testM_sparse.colptr;
rowind = testM_sparse.rowval;
nzvals = testM_sparse.nzval;
nnz = length(nzvals);


nnz = convert(Int32,nnz);
colptr = convert(Array{Int32,1},colptr);
rowind = convert(Array{Int32,1},rowind);

Lnnz = [0]
Lnnz = convert(Array{Int32,1},Lnnz);
permout=zeros(nnodes);
permout = convert(Array{Int32,1},permout);
V = ccall((:selinv2julia,
        "/home/grad/wgong/Downloads/newSelInv/SelInv/EXAMPLES/julia.so"), Ptr{Cdouble}, (Cint,Cint,Ptr{Cint},Ptr{Cint},Ptr{Cdouble},Ptr{Cint},Ptr{Cint}),
       nnodes,nnz,colptr,rowind, nzvals,Lnnz,permout);

v2 = pointer_to_array(V, (3*Lnnz[1],));
nnzInv = Lnnz[1];
iInv = v2[1:nnzInv];
jInv = v2[(nnzInv+1):(2*nnzInv)];
iInv =  convert(Array{Int32,1},iInv);
jInv = convert(Array{Int32,1},jInv);

for i=1:nnzInv
    iInv[i] = permout[iInv[i]];
    jInv[i] = permout[jInv[i]];
end

testM_inv = sparse(iInv,jInv,
	  v2[(2*nnzInv+1):(3*nnzInv)])
testM_inv = testM_inv + transpose(testM_inv)
for i=1:nnodes
    testM_inv[i,i] = testM_inv[i,i]/2
end

return full(testM_inv)
end 



inv(testM)
