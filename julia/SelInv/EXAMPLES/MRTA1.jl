cd("/home/grad/wgong/Downloads/newSelInv3/SelInv/EXAMPLES/")
function SelInv(testM_sp,B_sp)

  testM_sparse= testM_sp+B_sp*0.000000001

nnodes=size(testM_sp,1)
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



#cd("S:/work/distributedComputing/code")
#include("MRAfunctions.jl")
using Distances  # for distance computation (pairwise())

## sparse cholesky decomposition (returns A=LL' for sparse A)
function sparsechol(A)
  try
    return cholfact(A)
  catch e
    #println("caught error: $e")
    succ=false; di=size(A,1); attempt=1; B=0;
    while succ==false
      A += 10.0^(-10+attempt)*speye(Float64,di)
      try
        B=cholfact(A)
        succ=true
      end
      attempt += 1
    end
    println("$attempt attempts")
    return B
  end
end


## ceiling function (for easy switching between julia versions)
function ceilfun(x)
  iceil(x)  # ceil(Int64,x)
end


## simulate stationary data on 1-D grid (function DHSimulate in R package ltsa)
function DHsim(r)
  m = length(r)
  N = 2^ceilfun(log2(m-1)) #2^ceil(Int64,log2(m-1))
  acvf = [r; zeros(N-m+1)]
  append!(acvf,acvf[end-1:-1:2])
  g = real(fft(acvf))
  if (any(g .< 0.0))
    println("Davies-Harte nonnegativity condition not valid")
    return NaN
  else
    Z = complex(randn(N-1),randn(N-1))
    Z2 = 2 + sqrt(2) * randn(2)
    Z = [Z2[1]; Z; Z2[2]; conj(reverse(Z))]
    X = real(bfft(sqrt(g).*Z))/sqrt(2*N)
    z = X[1:m]/sqrt(2)
    return z
  end
end



## convenience functions for sparse matrices (will be implemented in Julia 0.4)
rowvals(S::SparseMatrixCSC) = S.rowval
nzrange(S::SparseMatrixCSC, col::Integer) = S.colptr[col]:(S.colptr[col+1]-1)

## determine nonzero elements based on distance  ## need faster algorithm for 1-D and 2-D
function distsparse(locs1,locs2,maxdist)
  rowinds=zeros(Int64,0); colinds=Int64[]; vals=Float64[];
  for i=1:size(locs1,1)
    for j=1:size(locs2,1)
      d=sqrt((locs1[i]-locs2[j])^2)
      if d<=maxdist
        append!(rowinds,[i])
        append!(colinds,[j])
        append!(vals,[d])
      end
    end
  end
  return {rowinds,colinds,vals}
end



## return matrix symbolizing the nz entries
function nzmatrix(locs1,locs2,maxdist)
  temp=distsparse(locs1,locs2,maxdist)
  return sparse(temp[1],temp[2],ones(Float64,size(temp[1],1)))
end

## selected inversion of sparse spd matrix --- replace by fortran code https://math.berkeley.edu/~linlin/software.html
#function selInv(A,B)
#  inv(full(A)).*B  # B is matrix of ones and zeros that specifies the sparsity structure
#end

## calculate marginal predictive variances  --- rewrite using selInv!!!
function calcvars(Bp,Ppchol)
  np=size(Bp,1)
  varvec=Array(Float64,np)
  for i=1:np
    varvec[i] = sum(Bp[i,:]'.*full(Ppchol\Bp[i,:]'))
  end
  return varvec
end



# define distance measure appropriate for the domain of interest
function distfun(locs1,locs2)
  pairwise(Euclidean(),locs1',locs2')
end
# covariance function
covFunBasic=function(x)
  return  .95*exp(-x) + 0.05*(abs(x).<1e-10)
end;
function co(locs1,locs2)
    return covFunBasic(distfun(locs1,locs2)/.05) # pairwise calculates dists between columns
end;
function kanter(t) # tapering function of Kanter (1997)
# valid in R^3, minimizes the curvature at the origin, only has 2 derivatives at the
#    origin (not useful for matern with nu>2), for more details, see gneiting (2002), eq. 22
  twopit=2*pi*t + 1.0e-300;  # add tiny number to ensure kanter(0)=1
  R=(1-t).*sin(twopit)./twopit+(1/pi)*(1-cos(twopit))./twopit;
  R[t.>1]=0; R[t.==0]=1;
  return R
end
function taper(distances,range)
  return kanter(distances/range)
end


# simulate data
M0=10; r0=10; J=2;  # M0=3; r0=2; J=3;
n=r0*J^M0;
sigeps=.05
locs=[(1:n)/n;]; # locs=[1; n-1:-1:2; n;]/n;
# locs=locs[[1:10; end-10:end]]; n=length(locs)
domain=[locs[1];locs[end];];
y= (n<10000 ? vec(chol(co(locs,locs))*randn(n,1)) : DHsim(vec(co([locs[1];],locs))) );
z=vec(y+randn(n,1)*sigeps);
predlocs=linspace(domain[1],domain[2],100)
pred=true


# using Gadfly; plot(x=locs,y=z)


# MRTA set-up
M=M0;
nzrowprior=4; # desired nz elements per row in prior
rangemultiplier=nzrowprior/r0/2;
range0=(domain[2]-domain[1])*rangemultiplier; #range is the d in paper
if M==1
  rs=[r0*M0; n;];
  ranges=[range0; range0/J^M0*M0];
else
  rs=r0*J.^[0:M;];
  ranges=range0./J.^[0:M;];
end
knots=Array(Any,M+1)
for m=0:(M-1)
  knotdist=(domain[2]-domain[1])/rs[m+1]
  knots[m+1] = linspace(domain[1]+knotdist/2,domain[2]-knotdist/2,rs[m+1])
  # knots[m+1] = domain[1]+[1:rs[m+1];]/(rs[m+1]+1)*(domain[2]-domain[1])
end
knots[M+1]=locs;



#################   calculate prior for nz elements    ###########

# pre-calculate sparse matrices for prior
nzmat=Array(Any,M+1)
for m=0:M
  nzmat[m+1]=Array(SparseMatrixCSC,m+1)
  for l=0:m
    nzmat[m+1][l+1]=nzmatrix(knots[m+1],knots[l+1],ranges[l+1])
  end
end
if pred
  nzmatpred=Array(SparseMatrixCSC,M+1)
  for l=0:M; nzmatpred[l+1]=nzmatrix(predlocs,knots[l+1],ranges[l+1]); end
end

# pre-calculate sparse matrices for SelInv
nzmat_SelInv=Array(Any,M)
for m=0:M-1
  nzmat_SelInv[m+1]=nzmatrix(knots[m+1],knots[m+1],(2+2/J)*ranges[m+1])
end


###############   the following section is very slow and needs to be re-written     ############

# create prior quantities  ---- version with SelInv
Plist=Array(SparseMatrixCSC,M+1);  #a collection of sparse matrices
Pchol=Array(Any,M+1)
K=Array(SparseMatrixCSC,M)
V=Array(Any,M+1)
for m=0:M
  V[m+1]=deepcopy(nzmat[m+1])
  for l=0:m
    Vrows = rowvals(V[m+1][l+1])
    Vvals = nonzeros(V[m+1][l+1])
    for j=1:size(V[m+1][l+1],2)
      i=nzrange(V[m+1][l+1],j)
      dists=distfun(knots[m+1][Vrows[i]],[knots[l+1][j]])
      Vvals[i]=co(knots[m+1][Vrows[i]],[knots[l+1][j]]).*taper(dists,ranges[1])
      for k=1:l
        Vvals[i] -= V[m+1][k][Vrows[i],:]*(K[k]*(V[l+1][k][j,:]')) # this line is the main bottleneck!!
        Vvals[i] .*= taper(dists,ranges[k+1])
      end
    end
  end
  Plist[m+1]=.5*V[m+1][m+1]+.5*V[m+1][m+1]'
  Pchol[m+1]=sparsechol(Plist[m+1])
  if m<M; K[m+1]=SelInv(Plist[m+1],nzmat_SelInv[m+1]); end #see defination of selInv
end

##  we need to port the SelInv Fortran code
## also, sparse matrix multiplication seems to incur substantial overhead, and so repeated calls to
##   sparse matrix multiplication need to be reduced as much as possible
## consider the following example:
test1=sparse(hcat(ones(1000,3),zeros(1000,100000)))
@time test2=test1';
function loopmult(test1,test2); for i=1:size(test1,1); test1[i,:]*test2; end; end;
function loopmult(test1,test2); for i=1:1000; test1[i,:]*test2; end; end;
@time test1*test2;
@time loopmult(test1,test2) # this does the same as the line above, but is 30 times slower!

########    end of "slow" section ###########


B=V[M+1][1]
for m=1:M
  B=hcat(B,V[M+1][m+1])
end
P=Plist[1] #sparse(Pchol[1])*sparse(Pchol[1])'
for m=1:M
  P=blkdiag(P,Plist[m+1]) #sparse(Pchol[m+1])*sparse(Pchol[m+1])')
end

if pred
  Vp=deepcopy(nzmatpred)
  for l=0:M
    Vrows = rowvals(Vp[l+1])
    Vvals = nonzeros(Vp[l+1])
    for j=1:size(Vp[l+1],2)
      i=nzrange(Vp[l+1],j)
      dists=distfun(predlocs[Vrows[i]],[knots[l+1][j]])
      Vvals[i]=co(predlocs[Vrows[i]],[knots[l+1][j]]).*taper(dists,ranges[1])
      for k=1:l
        Vvals[i] -= Vp[k][Vrows[i],:]*(K[k]*(V[l+1][k][j,:]'))
        Vvals[i] .*= taper(dists,ranges[k+1])
      end
    end
  end
  Bp=Vp[1];  for m=1:M; Bp=hcat(Bp,Vp[m+1]); end
end



#################    inference    ####################
Ppost = P + B'*B/sigeps^2
# Ppost=Ppost[end:-1:1,end:-1:1]; P=P[end:-1:1,end:-1:1]; B=B[:,end:-1:1]; # reordering
Ppchol=cholfact(Ppost)
Bpz = B'*z
KtBpz = Ppchol\Bpz


## likelihood
logdetPrior = sum([logdet(x) for x in Pchol])
logdetPost = logdet(Ppchol)
neg2loglik = logdetPost - logdetPrior + 2*n*log(sigeps) + sum(z.^2)/sigeps^2 - sum(KtBpz.*Bpz)/sigeps^4
loglik = -0.5*neg2loglik - n/2*log(2*pi)

## predictions
predmean = Bp*KtBpz/sigeps^2
predvar = calcvars(Bp,Ppchol)
predsd = sqrt(predvar)

using Gadfly
subset=1:100:n # 1:n
plot( layer(x=locs[subset],y=y[subset],Geom.point,Theme(default_color=color("red"))),
     layer(x=predlocs,y=predmean,Geom.line,Theme(default_color=color("blue"))),
     layer(x=predlocs,y=predsd,Geom.line,Theme(default_color=color("green"))))




#################   plotting and consideration of results    ###########
Pfact=sparsechol(P);
CM=B*(Pfact\(B'))

[nnz(P[:,i]) for i in 1:size(P,2)]
using Compose, Gadfly
subset=1:n^2;
dists=vec(pairwise(Euclidean(),locs',locs'))
plot(x=dists[subset],y=vec(CM)[subset])
spy(Bp)
spy(CM)
C0=co(locs,locs);
spy(C0)


nnz(B)/size(B,1)
nnz(Ppost)/size(Ppost,1)  # nz per row
nnz(Ppost)^2/size(Ppost,1)  # computational complexity
nnz(Ppost)/prod(size(Ppost))  # sparsity



#################   crude prior calculation (for comparison)    ###########
Plist=Array(Any,M+1)
Pchol=Array(Any,M+1)
V=Array(Any,M+1)
for m=0:M
  V[m+1]=Array(SparseMatrixCSC,m+1)
  for l=0:m
    V[m+1][l+1]=sparse(co(knots[m+1],knots[l+1]).*taper(knots[m+1],knots[l+1],ranges[1]))
    for k=1:l
      quadform=V[m+1][k]*(Pchol[k]\(V[l+1][k]'))
      V[m+1][l+1]=sparse((V[m+1][l+1]-quadform).*taper(knots[m+1],knots[l+1],ranges[k+1]))
    end
  end
  Plist[m+1]=.5*V[m+1][m+1]+.5*V[m+1][m+1]'
  Pchol[m+1]=sparsechol(Plist[m+1]) #,perm=1:size(V[m+1][m+1],1))
end
B=V[M+1][1]
for m=1:M
  B=hcat(B,V[M+1][m+1])
end
P=Plist[1] #sparse(Pchol[1])*sparse(Pchol[1])'
for m=1:M
  P=blkdiag(P,Plist[m+1]) #sparse(Pchol[m+1])*sparse(Pchol[m+1])')
end


########### simple sparse example  ####
A=sparse([1; 1; 2; 3; 3;],[1; 3; 2; 1; 3],[1:5;])
full(A)
rows = rowvals(A)
vals = nonzeros(A)
m, n = size(A)
for j = 1:n
  i=nzrange(A,j)
  row = rows[i]
  val = vals[i] # can modify column of A by modifying this
  # corresponding (full) indices are (row[i[1]],j),...,(row[i[end]],j)
      # perform sparse wizardry...
end
