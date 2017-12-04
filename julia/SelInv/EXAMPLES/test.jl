cd("/home/grad/wgong/Downloads/newSelInv3/SelInv/EXAMPLES/") ##not work

nnodes = 8;
nnz = 12;
testM = eye(nnodes)*4;
testM[2,1] = 1
testM[1,2] = 1
testM[3,2] = 1.2
testM[2,3] = 1.2

B = eye(nnodes)*8;
B[2,1] = 1
B[1,2] = 1
B[3,2] = 1
B[2,3] = 1
B[5,1]=1
B[6,2]=1
B[1,5]=1
B[2,6]=1

#testM=B*4
B_sp = sparse(B)
testM_sp= sparse(testM);

function SelInv(testM_sp,B_sp)
nnodes=size(testM_sp,1)
  testM_sparse= testM_sp+B_sp*0.000000001

colptr = testM_sparse.colptr;
rowind = testM_sparse.rowval;
nzvals = testM_sparse.nzval;
nnz = length(nzvals);

nnodes = convert(Int32,nnodes);
nnz = convert(Int32,nnz);
colptr = convert(Array{Int32,1},colptr);
rowind = convert(Array{Int32,1},rowind);

Lnnz = [0]
Lnnz = convert(Array{Int32,1},Lnnz);
nnzlplus = [0]
nnzlplus = convert(Array{Int32,1},nnzlplus);
permout=zeros(nnodes);
permout = convert(Array{Int32,1},permout);
DIAG=zeros(nnodes);
LDL_D = zeros(nnodes);
#/home/grad/wgong/Downloads/newSelInv3/SelInv/EXAMPLES/
V = ccall((:selinv2julia,
        "/home/grad/wgong/Downloads/newSelInv3/SelInv/EXAMPLES/julia.so"), Ptr{Cdouble}, (Cint,Cint,Ptr{Cint},Ptr{Cint},Ptr{Cdouble},
	Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cdouble},Ptr{Cdouble}),
       nnodes,nnz,colptr,rowind, nzvals,Lnnz,permout,nnzlplus,DIAG,LDL_D);

v2 = pointer_to_array(V, (3*Lnnz[1]+nnzlplus[1],));
nnzInv = Lnnz[1];
iInv = v2[1:nnzInv];
jInv = v2[(nnzInv+1):(2*nnzInv)];
iInv =  convert(Array{Int32,1},iInv);
jInv = convert(Array{Int32,1},jInv);
LDL_L =  v2[(2*nnzInv+1):(3*nnzInv)];
invElement =  v2[(3*nnzInv+1):(3*nnzInv+nnzlplus[1])];

inv_pos = sparse(iInv,jInv,1);

#DIAG = [1.0,2.0,5.0,10.0,15.0,20.0,21.0,22.0];

newLNZ = zeros(nnzInv);
for i=1:nnzInv
    if i == 1
       k = DIAG[jInv[i]];
    elseif jInv[i]>jInv[i-1]
       k = DIAG[jInv[i]];
    end
    newLNZ[i] = invElement[k];
    k = k + 1;
end


for i=1:nnzInv
    iInv[i] = permout[iInv[i]];
    jInv[i] = permout[jInv[i]];
end


## get the LDL component
LDL_L_matrix = sparse(iInv,jInv,LDL_L);
for i=1:nnodes
    LDL_L_matrix[i,i] = 1;
end
LDL_D_matrix = sparse([1:1:nnodes],[1:1:nnodes],LDL_D);

## reconstruct testM
testM_reconstruct = LDL_L_matrix*LDL_D_matrix*transpose(LDL_L_matrix);


testM_inv = sparse(iInv,jInv,newLNZ);
testM_inv = testM_inv + transpose(testM_inv);
for i=1:nnodes
    testM_inv[i,i] = testM_inv[i,i]/2
end

  return(testM_inv)
end



SelInv(testM_sp,B_sp)

full(testM_inv)
inv(testM)

LDL_D_matrix
LDL_L_matrix
