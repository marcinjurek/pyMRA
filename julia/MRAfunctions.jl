############################ function for estimation
## function to calculate MRAtaper likelihood
function MRAtaperll(knots,locs,z,r0,J,M,ranges,theta,sm,dmat,dmat1,dmat_SelInv,spbdS,spbdS1)#distmatrix,MRAprior_taper,MRApriorS,MRAloglik,
  sig=theta[1]
  prange=theta[2]
  sigeps=0.000000001+theta[3]
  print(theta)
  n=length(z)

  (B,V,Vt,P,K,LogDet)=MRAprior_taper(dmat,dmat1,dmat_SelInv,knots,ranges,sig,prange,sm,sigeps)#nzmat,nzmat1 # output K for prediction
     if knots[M+1]==locs
       loglik_MRTA=MRAloglikWE(P,B,LogDet,sigeps,z)
     else ##if Q_M not equal to locs then need
       BS=MRApriorS_taper(spbdS,spbdS1,locs,V,Vt,K,sig,prange,sm,sigeps,ranges)
       loglik_MRTA=MRAloglikWE(P,BS,LogDet,sigeps,z)
     end
    return loglik_MRTA
end
## function to calculate MRAblock likelihood
function MRAblockll(knots,locs,z,r0,J,M,theta,sm,spbd,spbd1,sp_Selinv,spbdS,spbdS1)#distmatrix,MRAprior_taper,MRApriorS,MRAloglik,
  sig=theta[1]
  prange=theta[2]
  sigeps=0.0000000001+theta[3]
  n=length(z)
  print(theta)
  (B,V,Vt,P,K,LogDet)=MRAprior_block(spbd,spbd1,sp_Selinv,sig,prange,sm,sigeps)
  if knots[M+1]==locs
    loglik_MRBA=MRAloglikWE(P,B,LogDet,sigeps,z)
  else ##if Q_M not equal to locs then need
    BS=MRApriorS(spbdS,spbdS1,locs,V,Vt,K,sig,prange,sm,sigeps)
    loglik_MRBA=MRAloglikWE(P,BS,LogDet,sigeps,z)
  end
  return loglik_MRBA
end

# loglikelihood for 1RA when using same set of knots from MRA
function OneRAtaperll(NewKnots,locs,z, r01,r0,J,M, Newranges,theta,sm,dmat,dmat1,dmat_SelInv)
  sig=theta[1]
  prange=theta[2]
  sigeps=0.000000001+theta[3]
  n=length(z)
  (B,V,Vt,P,K,LogDet)=MRAprior_taper(dmat,dmat1,dmat_SelInv,NewKnots,Newranges,sig,prange,sm,sigeps)#output K for prediction
    #if knots[M+1]==locs
      loglik_1RTA=MRAloglikWE(P,B,LogDet,sigeps,z)
    #else ##if Q_M not equal to locs then need
    #  BS=MRApriorS_taper(spbdS,spbdS1,locs,V,Vt,K,sig,prange,sm,sigeps,Newranges)
    #  loglik_1RTA=MRAloglikWE(P,BS,LogDet,sigeps,z)
    #end
    return loglik_1RTA
end
## loglikelihood for 1RA when use same set of knots from MRA
function OneRAblockll(NewKnots,locs,z, r01,r0,J,M, theta,sm,spbd,spbd1,sp_Selinv)

  sig=theta[1]
  prange=theta[2]
  sigeps=0.00001+theta[3]
  n=length(z)
  (B,V,Vt,P,K,LogDet)=MRAprior_block(spbd,spbd1,sp_Selinv,sig,prange,sm,sigeps)
  #if knots[M+1]==locs
    loglik_1RBA=MRAloglikWE(P,B,LogDet,sigeps,z)
  #else ##if Q_M not equal to locs then need
  #  BS=MRApriorS(spbdS,spbdS1,locs,V,Vt,K,sig,prange,sm,sigeps)
  #  loglik_1RBA=MRAloglikWE(P,BS,LogDet,sigeps,z)
  #end
  return  loglik_1RBA
end
## likelihood function for optimization, taper
function nloglikT(knots,locs,z,r0,J,M,ranges,logtheta,sm,dmat,dmat1,dmat_SelInv,spbdS,spbdS1)
    theta=exp(logtheta)
    nll=-MRAtaperll(knots,locs,z,r0,J,M,ranges,theta,sm,dmat,dmat1,dmat_SelInv,spbdS,spbdS1)#MRAtaperll(knots,locs,z,r0,J,M,ranges,alist,sm,)
  return nll
end
## likelihood function for optimization, block
function nloglikB(knots,locs,z,r0,J,M,logtheta,sm,spbd,spbd1,sp_Selinv,spbdS,spbdS1)
    theta=exp(logtheta)
    llb=-MRAblockll(knots,locs,z,r0,J,M,theta,sm,spbd,spbd1,sp_Selinv,spbdS,spbdS1)#MRAblockll(knots,locs,z,r0,J,M,theta,sm)
   return llb
end

## likelihood function for optimization, 1taper
function nloglik1RAT(NewKnots,locs,z, r01,r0,J,M,Newranges,logtheta,sm,dmat,dmat1,dmat_SelInv)
  theta=exp(logtheta)
  return -OneRAtaperll(NewKnots,locs,z, r01,r0,J,M,Newranges,theta,sm,dmat,dmat1,dmat_SelInv)
end
## likelihood function for optimization, 1block
function nloglik1RAB(NewKnots,locs,z, r01,r0,J,M, logtheta,sm,spbd,spbd1,sp_Selinv)
  theta=exp(logtheta)
  return -OneRAblockll(NewKnots,locs,z, r01,r0,J,M, theta,sm,spbd,spbd1,sp_Selinv)
end

#############################function for prediction
## sparsity structure matrices for prediction
function TaperMat_p(distmatrix,predlocs,knots,ranges)
  M=length(knots)-1
 spbdS_p=Array(Any,M+1)
 spbdS1_p=Array(Any,M+1)
   for m=0:M
     spbdS_p[m+1]=distmatrix(predlocs,knots[m+1],ranges[m+1])
     spbdS1_p[m+1]=distmatrix(knots[m+1],predlocs,ranges[m+1])
   end
   return (spbdS_p,spbdS1_p)
end
#=
## calculate marginal predictive variances  --- rewrite using selInv!!!
function calcvars(BSp,Ppchol)
  np=size(Bp,1)
  varvec=Array(Float64,np)
  for i=1:np
    varvec[i] = sum(BSp[i,:]'.*full(Ppchol\BSp[i,:]))
  end
  return varvec
end
=#


############################################# basic MRA functions
using Iterators #for product()
function D2knots(r0,J,M,domainX,domainY)
  n=r0*J^M
  knotsx=Array(Array,M+1)
  knotsy=Array(Array,M+1)
  knotsxy=Array(Array,M+1)
  rs=r0*J.^[0:M;];
  for m=0:M#(M-1)
    knotdistX=(domainX[2]-domainX[1])/rs[m+1]
    knotdistY=(domainY[2]-domainY[1])/rs[m+1]
    knotsx[m+1] = collect(linspace(domainX[1]+knotdistX/2,domainX[2]-knotdistX/2,rs[m+1]))
    knotsy[m+1] = collect(linspace(domainY[1]+knotdistY/2,domainY[2]-knotdistY/2,rs[m+1]))
      knotsxy[m+1]=collect(product(knotsx[m+1],knotsy[m+1]))
  end
#  knotsxy[M+1]=collect(product((0.5:1:n-0.5)/n,(0.5:1:n-0.5)/n))
  #Remove overlapped knots accross resolution
  knotsxy_u=Array(Array,M+1)
  knotsxy_all=Array(Array,M+1)
  knotsxy_all[1]=knotsxy[1]
  knotsxy_u[1]=knotsxy[1]
  for m=1:M
      knotsxy_all[m+1]=vcat(knotsxy_all[m],knotsxy[m+1])
      knotsxy_u[m+1]=setdiff(knotsxy[m+1],knotsxy_all[m])
  end
  #convert array of tuple to array of array
  knotsA=Array(Array,M+1)
  for m=0:M
  knotsA[m+1]=reinterpret(Float64, knotsxy_u[m+1], (2,length(knotsxy_u[m+1])))'
  end
  return knotsA
end

# Matern cov function with smoothness 0.5,1.5 or 2.5
covFunBasic=function(x,prange,sm) #sm= 0.5,1.5 or 2.5
  if sm==0.5
    cov=exp(-x/prange)
  else
    ind=ifelse(sm==2.5,1,0)
    cov=(1 + sqrt(2*sm)*x/prange+ ind*(2*sm*x.^2)/(3*prange.^2)).*exp(-sqrt(2*sm)*x/prange)
  end
  return  cov
end;
covFunWE=function(x,sig2,prange,sm,sigeps2)
        sig2*covFunBasic(x,prange,sm) + sigeps2*(abs(x/prange).==0)#<1e-40)#(1-sigeps)*
      end
##covariance function given two locations
function co(locs1,locs2,sig2,prange,sm,sigeps2)
    return covFunWE(distfun(locs1,locs2),sig2,prange,sm,sigeps2) # pairwise calculates dists between columns
end;

#=
function co(locs1,locs2,prange,sm,sigeps)
    return covFunWE(distfun(locs1,locs2),prange,sm,sigeps) # pairwise calculates dists between columns
end;
=#

## prior matrices function for MRAtaper
function MRApriorS_taper(spbdS,spbdS1,locs,knots,V,Vt,K,sig,prange,sm,sigeps,ranges)
  M=length(spbdS)-1
  VS=deepcopy(spbdS)#nzmat
  VS1=deepcopy(spbdS1)#nzmat
  A=Any
  for m=0:M
        dists = spbdS[m+1]
        VS[m+1] = co_spmat(dists,sig,prange,sm,sigeps).*taper_spmat(dists,ranges[1])
         for k=1:m
           A=K[k]*VS1[k]
            Vrows = (VS[m+1]).rowval #updated whenever Vals changed
            Vvals = (VS[m+1]).nzval #updated whenever Vals changed
            for j=1:size(VS[m+1],2)
              i=nzrange(VS[m+1],j) #return the range of indices of nz of a sp mat column
              #Vvals[i] -= VS[k][Vrows[i],:]*(K[k]*(V[m+1][k][j,:]))
              dists1=distfun(locs[Vrows[i],:],knots[m+1][j,:])
              Vvals[i] -= sum(broadcast(.*,A[:,Vrows[i]],Vt[m+1][k][:,j]),1)'
                Vvals[i] =   Vvals[i] .* taper(dists1,ranges[k+1])
            end
          end
          VS1[m+1]= VS[m+1]'
    end
#Stack Basis matrices and  precision matrices
    BS=VS[1]
    for m=1:M
      BS=hcat(BS,VS[m+1])
    end
   return BS
  end

##prior function for MRAblock (no tapering)
  function MRApriorS(spbdS,spbdS1,V,Vt,K,sig,prange,sm,sigeps)
    M=length(spbdS)-1
    VS=deepcopy(spbdS)#nzmat
    VS1=deepcopy(spbdS1)#nzmat
    A=Any
    for m=0:M
          dists = spbdS[m+1]
          VS[m+1] = co_spmat(dists,sig,prange,sm,sigeps)#.*taper_spmat(dists,ranges[1])
           for k=1:m
             A=K[k]*VS1[k]
              Vrows = (VS[m+1]).rowval #updated whenever Vals changed
              Vvals = (VS[m+1]).nzval #updated whenever Vals changed
              for j=1:size(VS[m+1],2)
                i=nzrange(VS[m+1],j) #return the range of indices of nz of a sp mat column
                #Vvals[i] -= VS[k][Vrows[i],:]*(K[k]*(V[m+1][k][j,:]))
                #dists1=distfun(locs[Vrows[i]],[knots[m+1][j]]')
                Vvals[i] -= sum(broadcast(.*,A[:,Vrows[i]],Vt[m+1][k][:,j]),1)'
                #  Vvals[i] =   Vvals[i] .* taper(dists1,ranges[k+1])
              end
            end
            VS1[m+1]= VS[m+1]'
      end
      #Stack Basis matrices and  precision matrices
      BS=VS[1]
      for m=1:M
        BS=hcat(BS,VS[m+1])
      end
     return BS
    end


##########################################################################
## Compute the logdet of B'B with input B
  function logdetB2(BS)
    LU=lufact(BS)
    logdetB2=2*sum(log(abs(diag(LU[:L])))) + 2*sum(log(abs(diag(LU[:U])))) - 2*sum(log(abs(LU[:Rs])))
    return logdetB2
  end

#=
  ## Compute the permuted cholesky factor of a sparse, symmetric positive definite matrix A.
  function cholfactPtL(A)
      n = size(A,1)
      F = cholfact(A)
      L0 = sparse(F[:L])
      is,js,es = findnz(L0)
      p = Base.SparseMatrix.CHOLMOD.get_perm(F)
      sparse(p[is], p[js], es, n, n)
  end
=#
## Compute the permuted cholesky factor of a sparse, symmetric positive definite matrix A.
function myPtL(F)
    L = sparse(F[:L])
    p = #Base.SparseArrays.CHOLMOD.get_perm(F)# SparseArrays for julia0.5
      Base.SparseMatrix.CHOLMOD.get_perm(F)#for julia 0.4
    return L[invperm(p),:]
end

################################################################################################
##loglik function when without measurement error
 function MRAloglik(P,BS,z)
   n=size(BS,1)
   Pchol = cholfact(P)
   Psqrt=myPtL(Pchol)
   z1=BS\z
   z2=Psqrt'*z1
   quad=sum(z2.^2)
   logD= logdetB2(BS) - logdet(Pchol)
   neg2loglik = quad + logD
   loglik = -0.5*neg2loglik - n/2*log(2*pi)
   return loglik
 end
##loglik function with full estimated cov matrix and mean
 function MRAloglikslow(P,BS,z)
      n=length(z)
      mu=0*ones(n)
      A=cholfact(P)\BS'
      Sigma_MRTA=BS*A
      Sigma_MRTA=0.5*Sigma_MRTA'+0.5Sigma_MRTA
      loglikMRTA=logpdf(MvNormal(mu,full(Sigma_MRTA)),z)
      return loglikMRTA
 end

###########################################################
#only works for model with measurement error
function MRAloglikWE(P,B,LogDet,sigeps,z)
  n=length(z)
  #posterior precison matrix
  Ppost = P + B'*B/sigeps #^2 # Ppost=Ppost[end:-1:1,end:-1:1]; P=P[end:-1:1,end:-1:1]; B=B[:,end:-1:1]; # reordering

  Ppchol=cholfact(Ppost)
  Bpz = B'*z
  KtBpz = Ppchol\Bpz

  ## log likelihood
  logdetPrior = sum(LogDet)
  #logdetPrior = 2*sum([logdet(full(x)) for x in Pchol])
  #logdetPrior1 = sum([logdet(x) for x in Pchol1])#because of magic cholfact output, don't need 2*
  logdetPost = logdet(Ppchol)
  neg2loglik = logdetPost - logdetPrior + n*log(sigeps) + sum(z.^2)/sigeps - sum(KtBpz.*Bpz)/sigeps^2# was 2*n*log(sigeps) + sum(z.^2)/sigeps^2 - sum(KtBpz.*Bpz)/sigeps^4
  loglik = -0.5*neg2loglik - n/2*log(2*pi)
  return loglik
end

 ################################################
 #loglik of the truth when n <=550000
 function DLloglik(dllik,locs,sigeps,z,prange,sm)
   n=length(z)
   if size(locs,2)==2
     r=sort(vec(co(locs[1,:]',locs,prange,sm,sigeps)),rev=true)
   else
     r=vec(co([locs[1];],locs,prange,sm,sigeps))#
   end
 r[1]=r[1]#+sigeps^2 #add measurement error on 0 distance covariance.but its already added in r above
 loglikDL= -0.5*dllik(r,z) - n/2*log(2*pi)
 return loglikDL
 end

#function to get New boundaries of each resolution(NewBound[m] to pick knots for m-th resolution)
 function Boundary(domain,J,M)
   numbd=Array(Int64,M+2) # number of bounds each resolution
   bound=Array(Any,M+2) # list of boundaries each resolution
   NewBound=Array(Vector,M+1)
   for m=0:M+1
       numbd[m+1]=J^m+1 # n + 1 boundaries give n subregions
       bound[m+1] = collect(linspace(domain[1],domain[2],numbd[m+1]))
       if m >0
         NewBound[m]= setdiff(bound[m+1],bound[m])
       end
   end
   return  NewBound
 end

#for regular spaced locs; assign the nearest loc on the right to each boundary as knot
#only works for r0=1, pick the first of R = J-1 knots.
#use bds from Boundary fct: all boundary values
function RightKnot(bds,locs,r0,M)
 kn=Array(Array,M+1)
   for m=0:M
     L=length(bds[m+1])
     kn[m+1]=Array(Float64,L)
     for l=1:L
     #find the nearest on the right
     kn[m+1][l]=locs[(locs .>= bds[m+1][l])][1]
     end
   end
   return kn
 end
#=
#plot boudaries and locs
using Gadfly
function PlotKnots(knots,locs)
 p=Gadfly.plot( #x=locs,y=ones(n),Geom.point,Coord.Cartesian(xmin=0,xmax=1) ,Theme(default_color=colorant"black")
 [ layer( x=knots[i], y=ones(length(knots[i]))*(i), Geom.point,Theme(default_color=colorant"red") ) for i in 1:M+1]... ,
   layer( x=locs,y=ones(n), Geom.point, Theme(default_color=colorant"black") ) )
 return p
end
=#
function Bounds(domain,J,M)
numbd=Array(Int64,M+1) # number of bounds each resolution
bds=Array(Any,M+1) # list of boundaries each resolution
for m=0:M
    numbd[m+1]=J^m + 1 # n+1 boundaries give n subregions
    bds[m+1] = collect(linspace(domain[1],domain[2],numbd[m+1]))
end
bounds=Array(Array,M+1) # bounds of each resolution
bounds[1]=domain'
for m=1:M
  bounds[m+1]=Array(Float64,numbd[m+1]-1,2)
  for l=0:numbd[m+1]-2
    bounds[m+1][l+1,:]=bds[m+1][(1:2).+l]'
  end
end
return bounds
end

function Bounds2D(domainX,domainY,J,M)
numbd=Array(Int64,M+1) # number of bounds each resolution
bdsX=Array(Any,M+1) # list of boundaries each resolution
bdsY=Array(Any,M+1) # list of boundaries each resolution
for m=0:M
    numbd[m+1]=J^m + 1 # n+1 boundaries give n subregions
    bdsX[m+1] = collect(linspace(domainX[1],domainX[2],numbd[m+1]))
    bdsY[m+1] = collect(linspace(domainY[1],domainY[2],numbd[m+1]))
end
boundsX=Array(Array,M+1) # bounds of each resolution
boundsY=Array(Array,M+1) # bounds of each resolution
boundsX[1]=domainX'
boundsY[1]=domainY'
for m=1:M
  boundsX[m+1]=Array(Float64,numbd[m+1]-1,2)
  boundsY[m+1]=Array(Float64,numbd[m+1]-1,2)
  for l=0:numbd[m+1]-2
    boundsX[m+1][l+1,:]=bdsX[m+1][(1:2).+l]'
    boundsY[m+1][l+1,:]=bdsY[m+1][(1:2).+l]'

  end
end
return (boundsX,boundsY)
end

## NewBounds is for 1RA block, #regions determined by #NewKnots/r0.
function NewBounds(domain,NewKnots,r0,M1)
  bounds=Array(Array,M1+1) # bounds of each resolution
  numbd=ones(Int64,2)
  #numbd[1]=Int64(round(length(NewKnots[1])/r0,0))
  numbd[2]=Int64(round(length(NewKnots[2])/r0,0))
  #rs=length.(NewKnots)
  #numbd=Int64.(round.(rs./r0,0))
  numbd[1]=2#Int64(round(numbd[2]/2,0))
  bds=Array(Any,M1+1) # list of boundaries each resolution
  for m=0:M1
      bds[m+1] = collect(linspace(domain[1],domain[2],numbd[m+1]+1)) # n+1 boundaries give n subregions
  end
  for m=0:M1
    bounds[m+1]=Array(Float64,numbd[m+1],2)
    for l=0:numbd[m+1]-1
      bounds[m+1][l+1,:]=bds[m+1][(1:2).+l]'
    end
  end
  return bounds
end

## generate knots for any J and R but when R!=J-1, generate overlapped knots across resolution
function Knots0(bounds,r0,J,M)
###calculate knots depend on bounds(J) and R###
kns =round(Int64,r0*J.^(0:M)) # number of knots each resolution
knots=Array(Any,M+1)
knotsall=Array(Float64,0)
for m=0:M
  knots[m+1]=Array(Float64,kns[m+1])
  for l=0:Int(kns[m+1]/r0-1)
    knots[m+1][l*r0+(1:r0)]=collect(linspace(bounds[m+1][l+1,1],bounds[m+1][l+1,2],2+r0))[2:r0+1]
  end
  knotsall = append!(knotsall,knots[m+1]) # cancatenate all knots into last resolution
end
  return (knots,knotsall)
end

### generate knots for any J and R=J-1(guarenteed no overlap), only take the first r0 knots of each subregion
#only works for r0 <= J-1, basicly pick a subset of R = J-1 knots.
#use bounds from Bounds function: pairwise boundaries
function Knots(bounds,r0,J,M)
  R=J-1
###calculate knots depend on bounds(J) and R###
kns =round(Int64,r0*J.^(0:M)) # number of knots each resolution
knots=Array(Any,M+1)
knotsall=Array(Float64,0)
for m=0:M
  knots[m+1]=Array(Float64,kns[m+1])
  for l=0:Int(kns[m+1]/r0-1)
    knots[m+1][l*r0+(1:r0)]=collect(linspace(bounds[m+1][l+1,1],bounds[m+1][l+1,2],2+R))[2:r0+1]
  end
  knotsall = append!(knotsall,knots[m+1]) # cancatenate all knots into last resolution
end
  return (knots,knotsall)
end

#Move Bounds slightly so that avoid overlapped knots across resolution
function Knots1(bounds,r0,J,M)
###calculate knots depend on bounds(J) and R###
kns =round(Int64,r0*J.^(0:M)) # number of knots each resolution
knots=Array(Any,M+1)
#overlap=Array(Vector,M)
knotsdist=Any
knots[1]=collect(linspace(bounds[1][1,1],bounds[1][1,2],r0+2))[2:r0+1]
knotsall=deepcopy(knots[1])
for m=1:M
  #m=3
  knots[m+1]=Array(Float64,kns[m+1])
  for l=0:size(bounds[m+1],1)-1
    knotsdist= (bounds[m+1][l+1,2]-bounds[m+1][l+1,1])/(r0+1)
    knots[m+1][l*r0+(1:r0)]=collect(linspace(bounds[m+1][l+1,1]+knotsdist/(2*m),bounds[m+1][l+1,2]+knotsdist/(2*m),2+r0))[2:r0+1]
  end
  #overlap[m]=intersect(knots[m+1],knotsall)
  knotsall = append!(knotsall,knots[m+1]) # cancatenate all knots into last resolution
end
  return (knots,knotsall) ##knotsall is before sorting
end

## uniformly generate knots for each region
function Knots2(bounds,r0,J,M)
###calculate knots depend on bounds(J) and R###
kns =round(Int64,r0*J.^(0:M)) # number of knots each resolution
knots=Array(Any,M+1)
#overlap=Array(Vector,M)
#knotsdist=Any
knots[1]=collect(linspace(bounds[1][1,1],bounds[1][1,2],r0+2))[2:r0+1]
knotsall=deepcopy(knots[1])
for m=1:M
  knots[m+1]=Array(Float64,kns[m+1])
  for l=0:size(bounds[m+1],1)-1
    #knotsdist= (bounds[m+1][l+1,2]-bounds[m+1][l+1,1])/(r0+1)
    knots[m+1][l*r0+(1:r0)]=rand(Uniform(bounds[m+1][l+1,1],bounds[m+1][l+1,2]),r0)
  end
  #overlap[m]=intersect(knots[m+1],knotsall)
  knotsall = append!(knotsall,knots[m+1]) # cancatenate all knots into last resolution
end
  return (knots,knotsall) ##knotsall is before sorting
end


function Knots3(domain,r0,J,M)
###calculate knots depend on bounds(J) and R###
kns =round(Int64,r0*J.^(0:M)) # number of knots each resolution
knotsdist=1./kns
knots=Array(Any,M+1)
bounds=Array(Float64,M+1,2)
bounds[:,1]=domain[1]-knotsdist/2
bounds[:,2]=domain[2]+knotsdist/2
knotsall=Array(Float64,0)
for m=0: M
  #bounds[m+1,1]=domain[1]-knotsdist[m+1]/(m+2)
  #bounds[m+1,2]=domain[2]+knotsdist[m+1]/(m+2) #middle knot overlap accross res.
  knots[m+1]=collect(linspace(bounds[m+1,1],bounds[m+1,2],kns[m+1]+2))[2:kns[m+1]+1]
  knotsall = append!(knotsall, knots[m+1])
end
return (knots,knotsall)
end



using Distributions
function GPll(co,locs,z,sig,prange,sm,sigeps)
  n=length(z)
  Sigma = co(locs,locs,sig,prange,sm,sigeps)
  #co(locs1,locs2,sig,prange,sm,sigeps)
  #C=co(locs,locs,prange,sm,0.05)
  mu=0*ones(n)
  loglikGP=logpdf(MvNormal(mu,Sigma),z)
  return loglikGP
end

## generate simulated data and knots together
function gkndat(r0,J,M,var,prange,sm,sigeps,seed)
  srand(seed)
  n=r0*J^M;
  locs=[(1:n)/n;]; # locs=[1; n-1:-1:2; n;]/n;
  domain=[locs[1];locs[end];];
  y= (n<10000 ? vec(chol(co(locs,locs,var,prange,sm,sigeps))'*randn(n,1)) : DHsim(vec(co([locs[1];],locs,var,prange,sm,sigeps))) );
  z=vec(y)

  ###Knots
  knots=Array(Any,M+1)
  rs=r0*J.^[0:M;];
  for m=0:(M-1)
    knotdist=(domain[2]-domain[1])/rs[m+1]
    knots[m+1] = collect(linspace(domain[1]+knotdist/2,domain[2]-knotdist/2,rs[m+1]))
    # knots[m+1] = domain[1]+[1:rs[m+1];]/(rs[m+1]+1)*(domain[2]-domain[1])
  end
  knots[M+1]=locs;

  return (knots,locs,z)
end
## function to generate simulated data and knots together in 2D domain
function gkndat2D(domainX,domaimY,r0,J,M,sm,sig,prange,sigeps,seed)
  n=r0*J^M
  knots=D2knots(r0,J,M,domainX,domainY)
  ##didn't find how to reinterpret by colum so have to transpose in the end
  #@time locs=reinterpret(Float64,collect(product((0.5:1:n-0.5)/n,(0.5:1:n-0.5)/n)),(2,n^2))' #was 0.5, 0.5
  grid=(0.5:1:n-0.5)/n
  locs=expandgrid(grid,grid)
  knots[M+1]=deepcopy(locs)

  #generate data with matern cov with m.e.
  #pre-compute the distacne matrix
  distM=distfun(locs,locs)
  y=  vec(chol(covFunWE(distM,sig,prange,sm,sigeps))'*randn(n^2,1))# DHsim not applicable for 2D, use RFsimulate in R
  #(n<10000 ? vec(chol(covFunWE(distM,prange,sm,sigeps))'*randn(n^2,1)): DHsim(vec(covFunWE(distM[:,1],prange,sm,sigeps))) );
  return (knots,locs,y)
end
# function to generate equi-distant data
function gdat(r0,J,M,sig2,prange,sm,sigeps,seed)
  srand(seed)
  n=r0*J^M;
  locs=[(1:n)/n;]; # locs=[1; n-1:-1:2; n;]/n;
  domain=[locs[1];locs[end];];
  y= (n<10000 ? vec(chol(co(locs,locs,sig2,prange,sm,sigeps))'*randn(n,1)) : DHsim(vec(co([locs[1];],locs,sig2,prange,sm,sigeps))) );
  z=vec(y)
  return (locs, z)
end

### function to generate knots given observed locs
function gkn(domain,r0,J,M,locs)
  ###Knots
  knots=Array(Any,M+1)
  rs=r0*J.^[0:M;];
  for m=0:(M-1)
    knotdist=(domain[2]-domain[1])/rs[m+1]
    knots[m+1] = collect(linspace(domain[1]+knotdist/2,domain[2]-knotdist/2,rs[m+1]))
    # knots[m+1] = domain[1]+[1:rs[m+1];]/(rs[m+1]+1)*(domain[2]-domain[1])
  end
  knots[M+1]=locs;

  return knots
end

##generate equal distant knots; put obs locs as the last res knots
function edknots(domain,r0,J,M,locs)
  knots=Array(Any,M+1)
  rs=r0*J.^[0:M;];
  for m=0:(M-1)
    knotdist=(domain[2]-domain[1])/rs[m+1]
    knots[m+1] = collect(linspace(domain[1]+knotdist/2,domain[2]-knotdist/2,rs[m+1]))
    # knots[m+1] = domain[1]+[1:rs[m+1];]/(rs[m+1]+1)*(domain[2]-domain[1])
  end
  knots[M+1]=locs;
  return knots
end

# generate for fixed n (i.e. fixed locs & knots)
# simulate data
function genknots(domain,r0,J,M,n,seed,sig2,sigeps,prange,sm)
  bs=Bounds(domain,J,M)
  knots=Array(Any,M+1)
  knotsall=Array(Float64,0)
  if (r0==J-1)
    (knots,knotsall)=Knots0(bs,r0,J,M)
  else
    (knots,knotsall)=Knots3(domain,r0,J,M) #Knots1 works for any J and r0
  end
  locs=sort(knotsall)# location equal to the union of all resolution knots
  srand(seed)
  #sigeps=0.0 #no measurement error (locs1,locs2,sig2,prange,sm,sigeps)
  y= ( n<10000 ? vec(chol(co(locs,locs,sig2,prange,sm,sigeps))'*randn(length(locs),1)) : DHsim(vec(co([locs[1];],locs,sig2,prange,sm,sigeps))) );
  z=vec(y)
  #Sigma = co(locs,locs)+ sigeps^2*eye(n) #with measurement error
  return (bs,knots,knotsall,locs,z) #need bs for MRAblock;#knotsall is before sorting used to checking overlap
end

## redistribute knots for 1RA based on knots of MRA
function oneRAknots(domain,knots,r0,M,n)
  NewKnots=Array(Any,2)
  NewKnots[1]=sort(collect(vcat(knots[1])))
  knotsunion=Array(Float64,0)
  if M>1
    NewKnots[2]=sort(collect([append!(knotsunion,knots[i]) for i=2:M+1 ])[1])
  else
    NewKnots[2]=sort(collect(knots[2]))
  end
  NewM=1
  NewJ=round(Int64,(n/r0-1))
  Newbs=Bounds(domain,NewJ,NewM)##regenerate bounds for 1RA, both r0 and J redefined.
  return (NewJ,NewM,NewKnots,Newbs)
end
#############################
### assign last resolution knots as the new 1RA knots: to sove the problem that B'z dimension not match problem when using oneRAknots
function oneRAknots1(domain,knots,r0,M,n)
  NewKnots=Array(Any,2)
  NewKnots[2]=sort(collect(vcat(knots[M+1])))
  knotsunion=Array(Float64,0)
  if M>1
    NewKnots[1]=sort(collect([append!(knotsunion,knots[i]) for i=1:M ])[1])
  else
    NewKnots[1]=sort(collect(knots[1]))
  end
  NewM=1
  NewJ=round(Int64,(n/r0))
  #Newbs=Bounds(domain,NewJ,NewM)##regenerate bounds for 1RA, both r0 and J redefined.
  return (NewJ,NewM,NewKnots)
end

#############################
# 1D version
### throw out the 1~M-1 resolution knots, assign the 0th and Mth res knots as NewKnots
function oneRAknots0(domain,knots,Newr0,r0)
  n=length(knots[end])
  NewKnots=Array(Any,2)
  if Newr0==r0
    NewKnots[1]=(sort(knots[1]))
  else
    d=(domain[2]-domain[1])/Newr0
    NewKnots[1]=collect(linspace(domain[1]+d/3,domain[2]-d/3,Newr0+2))[2:end-1]
  end
  NewKnots[2]=(sort(knots[end]))
  NewM=1
#Newr0=length(NewKnots[1])
  NewJ=ceil(Int64,(n/Newr0))
  return (NewJ,NewM,NewKnots)
end

#D2 version
### throw out the 1~M-1 resolution knots, assign the 0th and Mth res knots as NewKnots
function oneRAknots0D2(domain,knots,Newr0,r0)
  n=size(knots[end])[1]
  NewKnots=Array(Any,2)
  if Newr0==r0
    NewKnots[1]=deepcopy(knots[1])
  else
    d=(domain[2]-domain[1])/Newr0
    if size(knots[1],2)==2
      NewKnotsx=collect(linspace(domain[1]+d/3,domain[2]-d/3,Newr0+2))[2:end-1]
    #  NewKnotsy=collect(linspace(domain[1]+d/3,domain[2]-d/3,Newr0+2))[2:end-1]
    NewKnots[1]= expandgrid(NewKnotsx,NewKnotsx)
      #NewKnots1= collect(product(NewKnotsx,NewKnotsy))
      #NewKnots[1]=reinterpret(Float64, NewKnots1, (2,length(NewKnots1)))'
    else
      NewKnots[1]=collect(linspace(domain[1]+d/3,domain[2]-d/3,Newr0+2))[2:end-1]
    end
  end
  NewKnots[2]=deepcopy(knots[end])
  NewM=1
#Newr0=length(NewKnots[1])
  NewJ=ceil(Int64,sqrt(n)/Newr0)
  return (NewJ,NewM,NewKnots)
end

# redistribute MRA knots for 1RA in 2D domain
function oneRAknots2D(domainX,domainY,knots,Newr0,r0)
  n=size(knots[end])[1]
  NewKnots=Array(Any,2)
  if Newr0==r0
    NewKnots[1]=deepcopy(knots[1])
  else
    dx=(domainX[2]-domainX[1])/Newr0
    dy=(domainY[2]-domainY[1])/Newr0
    NewKnotsx=collect(linspace(domainX[1]+dx/3,domainX[2]-dx/3,Newr0+2))[2:end-1]
    NewKnotsy=collect(linspace(domainY[1]+dy/3,domainY[2]-dy/3,Newr0+2))[2:end-1]
    NewKnots[1]= expandgrid(NewKnotsx,NewKnotsy)
      #NewKnots1= collect(product(NewKnotsx,NewKnotsy))
      #NewKnots[1]=reinterpret(Float64, NewKnots1, (2,length(NewKnots1)))'
  end
  NewKnots[2]=deepcopy(knots[end])
  NewM=1
#Newr0=length(NewKnots[1])
  NewJ=ceil(Int64,sqrt(n)/Newr0)
  return (NewJ,NewM,NewKnots)
end

# redistribute MRA knots for 1RA
function oneRAknots0more(knots)
  n=length(knots[end])
  NewKnots=Array(Any,2)
  NewKnots[1]=sort(vcat(knots[1],knots[2]))
  NewKnots[2]=sort(knots[M+1])
  NewM=1
  Newr0=length(NewKnots[1])
  NewJ=round(Int64,(n/Newr0))
  return (Newr0,NewJ,NewM,NewKnots)
end

##Plot knots of 1RA
function plotKnots1RA(r0,J1,M1,NewKnots)
   ind=Array(UnitRange{Int64},2)
   l1=length(NewKnots[1])
   l2=length(NewKnots[2])
   ind[1]=1:l1
   ind[2]=l1+1:l2
   rs=[l1;l2]
   D=vcat([ createDFkn(m,ind,rs,NewKnots) for m in 0:M1 ]...)
   p = plot( D,
       layer(x=:x, y=:y, color=:m,  Geom.point) ,Coord.Cartesian(ymin=0,xmin=domain[1],xmax=domain[2])   )
return p
end

# Compute total NNZ of all resolutions
function TotNNZ(dmat)
  M=size(dmat,1)-1
  NNZtaper10=zeros(Int64,M+1,M+1)
  for m=0:M
    for l=0:m
      NNZtaper10[m+1,l+1]=nnz(dmat[m+1][l+1])
    end
  end
  return sum(NNZtaper10)
end

# Compute total NNZ of all resolutions
function TotNNZ1(dmat_SelInv)
  M=size(dmat_SelInv,1)-1
  NNZtaper10=zeros(Int64,M+1)
  for m=0:M
      NNZtaper10[m+1]=nnz(dmat_SelInv[m+1])
  end
  return sum(NNZtaper10)
end

# Compute total NNZ of all resolutions
function TotNNZmm(dmat)
  M=size(dmat,1)-1
  NNZtaper10=zeros(Int64,M+1)
  for m=0:M
    #for l=0:m
      NNZtaper10[m+1]=nnz(dmat[m+1][m+1])
    #end
  end
  return sum(NNZtaper10)
end

#if (VERSION >= v"0.5-")
###function to plot basis
using DataFrames,Gadfly
using Cairo,Fontconfig #to save pdf

function BasisCount(r0,J,M)
  rs=r0*J.^(0:(M))
  csumrs=cumsum(rs)
  ind=Array(UnitRange{Int64},M+1)
  ind[1]=1:csumrs[1]
  for m=1:M
    ind[m+1]=csumrs[m]+1:csumrs[m+1]
  end
  return(ind,rs)
end

function createDF(m::Int,ind,rs,BSf,locs)
    n=length(locs)
    i = ind[m+1]
    iln = rs[m+1]
    DataFrame(x=repmat(locs, iln), y=vec(BSf[:,i])+2*m, m="$m", i=repeat(i, inner=[n]),Taper="YES")
end

function createDF0(m::Int,ind,rs,BSf,locs)
    n=length(locs)
    i = ind[m+1]
    iln = rs[m+1]
    DataFrame(x=repmat(locs, iln), y=vec(BSf[:,i])+2*m, m="$m", i=repeat(i, inner=[n]),Taper="No")
end

function plotBasis(r0,J,M,BS)
  (ind,rs)=BasisCount(r0,J,M)
  BSf=full(BS) #need to convert BS to full matrix to make createDF work.
  D = vcat([ createDF(m,ind,rs,BSf,locs) for m in 0:M ]...)
  p = plot( D,
      layer(x=:x, y=:y, color=:m, group=:i, Geom.line),
      yintercept=[0,1,2,3], Geom.hline(color=colorant"darkgray", size=0pt),Guide.title("nzrp=$nzrowprior;r0=$r0;J=$J;M=$M")
      )
  return p
end
#create DataFrame for knots when plotting/testing non-overlapped knots generating
function createDFkn(m::Int,ind,rs,kn)
    i = ind[m+1]
    iln = rs[m+1]
    DataFrame(x=kn[m+1], y=ones(iln)*m, m="$m")#
end
#
function PlotKnots(doamin,r0,J,M,kn)
  (ind,rs)=BasisCount(r0,J,M)
  D=vcat([ createDFkn(m,ind,rs,kn) for m in 0:M ]...)
  p = plot( D,
      layer(x=:x, y=:y, color=:m,  Geom.point) ,Coord.Cartesian(ymin=0,xmin=domain[1],xmax=domain[2])   )
  return p
end

#end


##check overlapping  knots accross resolution
#the overlapped knots of each resolution with all the rest knots
function checkoverlap(r0,J,M,kn,knotsall) #knotsall is before sorting
  ind=BasisCount(r0,J,M)[1]
  overlap1=Array(Vector,M+1)
  for m=0:M
    catknots=knotsall[setdiff(1:end,ind[m+1])]  #splice!(locs,ind[m+1]) would modify the locs
    overlap1[m+1]=intersect(catknots,kn[m+1])
  end
 return overlap1
end


##To find effective range of exp cov
function effrange(sig,prange,sm)
  ln=1000
  dvec = collect(1:ln)/ln
  covec=covFunWE(dvec,sig,prange,sm,0)
  return dvec[find(covec.>= 0.05)][end]
end


##find the minimum distance between locs (sorted locs)
function mdist(locs)
ln=size(locs,1)
d=zeros(Float64,ln-1)
for i=1:ln-1
  d[i]=locs[i+1]-locs[i]
end
return sort(d)
end



function OneRAtaper(oneRAknots0D2,domain, locs,z,sig,prange, sm,sigeps, knots, r01, r0, ranges, rate,k)
  (J1,M1,NewKnots)=oneRAknots0D2(domain,knots,r01,r0)# r01=r0 indicate r0 the same
  if (size(locs,2)==1)
  Newranges=[ranges[1],ranges[end]].*rate ###this one not work for 2D
  J2=J1
  else
  J2=J1^2
  Newranges=TaperRangesNew(sig,prange,sm,r0,J2,M1,k)
  end
  (loglik1RTA,loglik_1RTA,T1RTA,TNNZT1)= MRAtaperWE(Newranges,NewKnots,locs,z,r01,J2,M1,sig,prange,sm,sigeps)
return (loglik_1RTA,T1RTA,TNNZT1)#,Newranges[1],Newranges[2])
end

function OneRAblock(oneRAknots0D2,domain, locs,z,sig,prange, sm, sigeps, knots, r01, r0)
    (J1,M1,NewKnots)=oneRAknots0D2(domain,knots,r01,r0)# 0 indicate r0 the same
    Newbs=Bounds(domain,J1,M1)
      if size(locs,2)==1
        J2= J1
    else
      J2=J1^2
    end
    (loglik_1RBA,T1RBA,TNNZB1)= MRAblockWE(Newbs,NewKnots,locs,z,r01,J2,M1,sig,prange,sm,sigeps)
    return  (loglik_1RBA,T1RBA,TNNZB1)
  end


  function createDF(m::Int,knots)
    DataFrame(x=knots[m+1][:,1], y=knots[m+1][:,2], m="$m")
  end
  function plotknots2D(knots)
    M=length(knots)-1
    DF=vcat([ createDF(m,knots) for m in 0:M]...)
    p = plot( DF,layer(x=:x, y=:y, color=:m,  Geom.point),Theme(default_point_size=2pt))
         #,Coord.Cartesian(ymin=0,xmin=domain[1],xmax=domain[2])
    return p
  end
function plotknots2D1res(knots)
  M=0#length(NewKnots[1])-1
  DF=vcat([ createDF(m,NewKnots) for m in 0:M]...)
  p = plot( DF,
      layer(x=:x, y=:y, color=:m,  Geom.point) #,Coord.Cartesian(ymin=0,xmin=domain[1],xmax=domain[2])
      )
  return p
end


  ### returns 2-column matrix containing all combinations of two vectors
  function expandgrid(vec1,vec2)
    combs=Array(typeof(vec1[1]),length(vec1)*length(vec2),2);
    counter=1;
    for x=1:length(vec1)
      for y=1:length(vec2)
        combs[counter,:]=[vec1[x];vec2[y];];
        counter += 1;
      end
    end
    return combs
  end

##################################################################################################################
# prediction for block version
function pred_block(z,predlocs,theta_hatB,knots,sm,boundsX,boundsY,spbd,spbd1,sp_Selinv,spbdS,spbdS1,spbdS_p,spbdS1_p)
  #theta_hatT=[0.95,0.05,0.05]
  sig2=theta_hatB[1]
  prange=theta_hatB[2]
  sigeps2=theta_hatB[3]
  ##predict
  (B,V,Vt,P,K,LogDet)= MRAprior_block(spbd,spbd1,sp_Selinv,sig2,prange,sm,sigeps2);
  #(spbdS_p,spbdS1_p)=BlockMat_p(spBd,knots,predlocs,boundsX,boundsY)
  ####if Q_M not equal to locs then need BS
  #BS = MRApriorS(spbdS,spbdS1,V,Vt,K,sig2,prange,sm,sigeps2) # because S=UQ_m, BS is the basis used for loglik and prediction
  BS=B # to avoid changing Ppost & Bpz below
  BSp = MRApriorS(spbdS_p,spbdS1_p,V,Vt,K,sig2,prange,sm,sigeps2)

  if (sigeps2==0)
    invBSz=BS\z # same as inv(full(BS))*z
    predmean=BSp*invBSz
  else
    Ppost = P +   BS'*  BS/sigeps2#^2 ## doesn't apply to sigeps=0
    Ppchol=cholfact(Ppost)
    Bpz =   BS'*z
    KtBpz = Ppchol\Bpz
    predmean = BSp*KtBpz/sigeps2#^2
  end

#marginal prediction err
  A=Ppchol\BSp'
  varM=BSp*A
  #predvar = calcvars(BSp,Ppchol) # not work, didn't figure out why
  predsd = sqrt(diag(varM))

  predresult=hcat(predlocs,predmean,predsd)
  return predresult
end
#find the bounds for x axis based on bounds(divide on both x and y), when J=2 first divide into 2 on x axis
function boundx(bounds)
  M=length(bounds)-1
  bdx=Array(Any,M+1) # list of boundaries each resolution
  bdx[1]=bounds[1]
  for m=2:M+1
    if iseven(m)==true
      bdx[m] = bounds[floor(Int32,m/2)+1]
    else
      bdx[m] = bdx[m-1]
    end
  end
  return bdx
end

#find the bounds for y axis based on bounds(divide on both x and y), when J=2 first divide into 2 on x axis
function boundy(bounds)
  M=length(bounds)-1
  bdy=Array(Any,M+1) # list of boundaries each resolution
  for m=1:M+1
    if isodd(m)==true
      bdy[m] = bounds[floor(Int32,m/2)+1]
    else
      bdy[m] = bdy[m-1]
    end
  end
  return bdy
end

#reorder the knots by blocks
function reod(knots,bdx,bdy)
  M=length(knots)-1
#knots[1]: 0 res
  m=0
  knots[m+1]=sortrows(knots[m+1],by=x->(x[1],x[2]))
#knots[2] and up : 1 res and up
  for m=1:M
    tmp=hcat(Float64[],Float64[])
    for i=1:size(bdx[m+1])[1] #or bdy, the same length
      for j=1:size(bdy[m+1])[1]
        #find knots[i] indexes within corresponding bdx[i] and bdy[i]
      ind=findin(knots[m+1],bdx[m+1][i,:],bdy[m+1][j,:])
      #Yind=find(bdy[k][i,1] .<= knots[k][:,2] .<= bdy[k][i,2])
      #order each by y axis and cancatenate them
      tmp=vcat(tmp,sortrows(knots[m+1][ind,:],by=x->(x[1],x[2])))
    end
    end
    knots[m+1]=tmp
  end
  return knots
end


function spB(knots1,knots2,bound)
  IJ=Array(Int64,0,2)
  dvals=ones(Float64,0)#empty vector
  for b=1:size(bound)[1]
    Is=find(bound[b,:][1] .< knots1 .<= bound[b,:][2])
    Js=find(bound[b,:][1] .< knots2 .<= bound[b,:][2])
    #D= distmat_d(knots1[Is],knots2[Js])
    #dvals_current=vec(D)
    #dvals_current[find(dvals_current.==0)]=e0
    #dvals=vcat(dvals,dvals_current)
    IJ_current = hcat(repeat(Is,inner=[1],outer=[length(Js)]),repeat(Js,inner=[length(Is)],outer=[1]))
    IJ=vcat(IJ,IJ_current)
  end
  return sparse(IJ[:,1],IJ[:,2],ones(Int64,size(IJ)[1]))
end

#=
function spBd(knots1,knots2,bound)
  e=1.0e-50
  m=size(knots1,1)
  n=size(knots2,1)
  IJ=Array(Int64,0,2)
  dvals=ones(Float64,0)#empty vector
  for b=1:size(bound)[1]
    #b=1
    Is=find(bound[b,:][1] .<= knots1 .<= bound[b,:][2])
    Js=find(bound[b,:][1] .<= knots2 .<= bound[b,:][2])
    D= distmat_d(knots1[Is],knots2[Js])
    dvals_current=vec(D)
    dvals_current[find(dvals_current.==0)]= e#[1-1] #make it exact 0 in julia 0.4.5 & 0.5.0
    dvals=vcat(dvals,dvals_current)
    IJ_current = hcat(repeat(Is,inner=[1],outer=[length(Js)]),repeat(Js,inner=[length(Is)],outer=[1]))
    IJ=vcat(IJ,IJ_current)
  end
  A=sparse(IJ[:,1],IJ[:,2],dvals,m,n)
  A.nzval[(A.nzval.==e)]=0
  return A
end
=#

#function to find which elements(2D) fall in the bound; reuturn the index of the elements
# need x in boundx and y in boundy
function findin(knots1,bdx,bdy)#knots1 is a list of 2D locations; bound is common for both x-axis and y-axis
  #bd=bound[1,:]
  Xind=find(bdx[1] .<= knots1[:,1] .<= bdx[2])
  Yind=find(bdy[1] .<= knots1[:,2] .<= bdy[2])
  return intersect(Xind,Yind)
end
#findin(knots1,bound[1,:])

function spBd(knots1,knots2,bound,boundy)#bound is for boundx, to avoid changing 1D case, use bound as name
  e=1.0e-50
  m=size(knots1,1)
  n=size(knots2,1)
  IJ=Array(Int64,0,2)
  dvals=ones(Float64,0) #empty vector
  if size(knots1,2)==2  #calc dist matrix for 2D locs
    for bx=1:size(bound)[1]
      for by=1:size(boundy)[1]
      Is=findin(knots1,bound[bx,:],boundy[by,:]) # indexes of knots1 that fall in the bound
      Js=findin(knots2,bound[bx,:],boundy[by,:]) # indexes of knots2 that fall in the boun
      #Is=find(bound[b,:][1] .<= knots1 .<= bound[b,:][2])
      #Js=find(bound[b,:][1] .<= knots2 .<= bound[b,:][2])
      D= distfun(knots1[Is,:],knots2[Js,:])#distmat_d(knots1[Is,:],knots2[Js,:])
      dvals_current=vec(D)
      dvals_current[find(dvals_current.==0)]= e#[1-1] #make it exact 0 in julia 0.4.5 & 0.5.0
      dvals=vcat(dvals,dvals_current)
      IJ_current = hcat(repeat(Is,inner=[1],outer=[length(Js)]),repeat(Js,inner=[length(Is)],outer=[1]))
      IJ=vcat(IJ,IJ_current)
    end
    end
  else
    for b=1:size(bound)[1]
      Is=find(bound[b,:][1] .<= knots1 .<= bound[b,:][2])
      Js=find(bound[b,:][1] .<= knots2 .<= bound[b,:][2])
      D= distmat_d(knots1[Is],knots2[Js])
      dvals_current=vec(D)
      dvals_current[find(dvals_current.==0)]= e#[1-1] #make it exact 0 in julia 0.4.5 & 0.5.0
      dvals=vcat(dvals,dvals_current)
      IJ_current = hcat(repeat(Is,inner=[1],outer=[length(Js)]),repeat(Js,inner=[length(Is)],outer=[1]))
      IJ=vcat(IJ,IJ_current)
    end
  end
  #full(sparse(IJ_current[:,1],IJ_current[:,2],dvals_current))
  A=sparse(IJ[:,1],IJ[:,2],dvals,m,n)
  A.nzval[(A.nzval.==e)]=0
  return A
end
## prior function for MRA block
function MRAprior_block(spbd,spbd1,sp_Selinv,sig,prange,sm,sigeps)
  M=length(spbd)-1
  Plist=Array(SparseMatrixCSC,M+1);  #a collection of sparse matrices
  LogDet=Array(Float64,M+1)
  K=Array(SparseMatrixCSC,M+1)
  V=deepcopy(spbd)#nzmat #
  A=Any
  Vt=deepcopy(spbd1)#nzmat1 #
  for m=0:M
     for l=0:m
        dists=spbd[m+1][l+1]
        V[m+1][l+1] = co_spmat(dists,sig,prange,sm,sigeps)#calc the value in V for k,j loop
        #m=3;l=1;k=1
       for k=1:l
          A=(K[k]*Vt[m+1][k]) ######
          Vrows = (V[m+1][l+1]).rowval #updated whenever Vals changed
          Vvals = (V[m+1][l+1]).nzval #updated whenever Vals changed
          for j=1:size(V[m+1][l+1],2)
            #j=1
            i=nzrange(V[m+1][l+1],j) #return the range of indices of nz of a sp mat column
            Vvals[i] -= sum(broadcast(.*,A[:,Vrows[i]],Vt[l+1][k][:,j]),1)'
            #Vvals[i] .*= spBv(knots[m+1][Vrows[i]],knots[l+1][j],bs[k])
          end
        end
      Vt[m+1][l+1] = V[m+1][l+1]'
    end
    Plist[m+1]=.5*V[m+1][m+1]+.5*V[m+1][m+1]';
    #K[m+1]=inv(full(Plist[m+1]))
    (K[m+1],LogDet[m+1])= SIplusD(Plist[m+1],sp_Selinv[m+1]);
  end

  #Stack Basis matrices and  precision matrices
  B=V[M+1][1]
  for m=1:M
    B=hcat(B,V[M+1][m+1])
  end
  P=Plist[1] #sparse(Pchol[1])*sparse(Pchol[1])'
  for m=1:M
    P=blkdiag(P,Plist[m+1]) #sparse(Pchol[m+1])*sparse(Pchol[m+1])')
  end

  return (B,V,Vt,P,K,LogDet)
end

## sparsity structure needed for prediction
function BlockMat_p(spBd,knots,locs,bs,bsy)
  M=length(knots)-1
spbdS_p=Array(Any,M+1)
spbdS1_p=Array(Any,M+1)
  for m=0:M
    spbdS_p[m+1]=spBd(locs,knots[m+1],bs[m+1],bsy[m+1])
    spbdS1_p[m+1]=spBd(knots[m+1],locs,bs[m+1],bsy[m+1])
  end
  return (spbdS_p,spbdS1_p)
end
## sparsity structure needed for prior calculation
function BlockMat(spBd,knots,locs,bs,bsy)
    M=length(bs)-1
    #spb=Array(Any,M+1)
    spbd=Array(Any,M+1)
    spbd1=Array(Any,M+1)
    spbdS=Array(Any,M+1)
    spbdS1=Array(Any,M+1)
    #if sigeps==0
        for m=0:M
        spbdS[m+1]=spBd(locs,knots[m+1],bs[m+1],bsy[m+1])
        spbdS1[m+1]=spBd(knots[m+1],locs,bs[m+1],bsy[m+1])
      end
  #  end
    for m=0:M
      #spb[m+1]=Array(SparseMatrixCSC,m+1)
      spbd[m+1]=Array(SparseMatrixCSC,m+1)
      spbd1[m+1]=Array(SparseMatrixCSC,m+1)

      #spb[m+1][1]= spB(knots[m+1],knots[1],bs[1])
      spbd[m+1][1]= spBd(knots[m+1],knots[1],bs[1],bsy[1])
      spbd1[m+1][1]= spBd(knots[1],knots[m+1],bs[1],bsy[1])
      #spbdS[m+1]=spBd(locs,knots[m+1],bs[m+1])
      #spbdS1[m+1]=spBd(knots[m+1],locs,bs[m+1])

      for l=1:m
        #spb[m+1][l+1]= spB(knots[m+1],knots[l+1],bs[l+1])
        spbd[m+1][l+1]= spBd(knots[m+1],knots[l+1],bs[l+1],bsy[l+1])
        spbd1[m+1][l+1]=spBd(knots[l+1],knots[m+1],bs[l+1],bsy[l+1])
      end
    end

    sp_Selinv=Array(Any,M+1)
    for m=0:M
      sp_Selinv[m+1]=spbd[m+1][m+1]#spB(knots[m+1],knots[m+1],bs[m+1])
    end
    return (spbd,spbd1,sp_Selinv,spbdS,spbdS1)
end

#=
function MRAblock(bs,knots,locs,z,r0,J,M,prange,sm)
    ## MRTA set-up
    #compute all related blcck matrices for MRAprior_block and MRApriorS
    (spbd,spbd1,spbdS,spbdS1,sp_Selinv)=BlockMat(spBd,knots,locs,bs)
    TNNZ=TotNNZ1(sp_Selinv)
  tic()
    #compute prior
    (B,V,Vt,P,K,LogDet)=MRAprior_block(spbd,spbd1,sp_Selinv,prange,sm)
    BS= MRApriorS(spbdS,spbdS1,V,Vt,K,prange,sm)
    T=toq()
    #compute logliklihood of MRAtaper
    loglik_MRBA=MRAloglik(P,BS,z)

    #compute the true GP logliklihood for comparison
    #loglikGP=GPll(co,locs,n,z)
    return (loglik_MRBA,T,TNNZ)
  end
=#
#loglikelihood of MRAblock with measurement error
  function MRAblockWE(bs,knots,locs,z,r0,J,M,sig,prange,sm,sigeps)
      ## MRTA set-up
       #bs=Bounds(domain,J,M)
      #compute all related blcck matrices for MRAprior_block and MRApriorS
      (spbd,spbd1,sp_Selinv)=BlockMat(spBd,knots,locs,bs,sigeps)[1:3]
      TNNZ=TotNNZ1(sp_Selinv)
      tic()
      #compute prior
      (B,V,Vt,P,K,LogDet)=MRAprior_block(spbd,spbd1,sp_Selinv,sig,prange,sm,sigeps)
      TB=toq()
      #compute logliklihood of MRAtaper
      loglik_MRBA=MRAloglikWE(P,B,LogDet,sigeps,z)

      #compute the true GP logliklihood for comparison
      #loglikGP=GPll(co,locs,n,z)
      return (loglik_MRBA,TB,TNNZ)
    end
########################################################################################3####################
#prediction of MRA taper
function pred_taper(distmatrix,z,predlocs,theta_hatT,knots,ranges,sm,dmat,dmat1,dmat_SelInv,spbdS,spbdS1,spbdS_p,spbdS1_p)
  sig=theta_hatT[1]
  prange=theta_hatT[2]
  sigeps=theta_hatT[3]
##predict
  (B,V,Vt,P,K,LogDet)=MRAprior_taper(dmat,dmat1,dmat_SelInv,knots,ranges,sig,prange,sm,sigeps)
  #(spbdS_p,spbdS1_p)= TaperMat_p(distmatrix,predlocs,knots,ranges)
  BSp=MRApriorS_taper(spbdS_p,spbdS1_p,predlocs,knots,V,Vt,K,sig,prange,sm,sigeps,ranges)
  Ppost = P +   B'*  B/sigeps#^2
  Ppchol=cholfact(Ppost)
  Bpz =   B'*z
  KtBpz = Ppchol\Bpz
  predmean = BSp*KtBpz/sigeps#^2
#marginal prediction err
  A=Ppchol\BSp'
  varM=BSp*A
  #predvar = calcvars(BSp,Ppchol) # not work, didn't figure out why
  predsd = sqrt(diag(varM))

  predresult=hcat(predlocs,predmean,predsd)
return predresult
end
#cd("/home/grad/wgong/Downloads/newSelInv3/SelInv/EXAMPLES/")
#A function to get selinv of matrix corresponding to certain sparsity spB
function SIplusD(spA,spB)
  spM=spA+1.0e-30*spB
  return SelInvD(spM)
end


# Selected inversion of sparse symmetrix matrix corresponding to nonzero entries.
# Output the selected inverse matrix and logdet of the input matrix
function SelInvD(testM_sparse)
#print(isposdef(testM_sparse))
nnodes=size(testM_sparse,1)
colptr = testM_sparse.colptr;
rowind = testM_sparse.rowval;
nzvals = testM_sparse.nzval;
nnz = length(nzvals);

nnodes = convert(Int32,nnodes);
nnz = convert(Int32,nnz);
colptr = convert(Array{Int32,1},colptr);
rowind = convert(Array{Int32,1},rowind);

Lnnz = Int32[0]
nnzlplus = Int32[0]
permout=zeros(Int32,nnodes);
#permout = convert(Array{Int32,1},permout);
DIAG=zeros(nnodes);
LDL_D = zeros(nnodes);
#=
V = ccall((:selinv2julia,"/home/grad/wgong/Spatial/newSelInv3/SelInvCopy/EXAMPLES/julia.so"),
        Ptr{Cdouble}, (Cint,Cint,Ptr{Cint},Ptr{Cint},Ptr{Cdouble},Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cdouble},Ptr{Cdouble}),
                        nnodes,nnz,colptr,rowind, nzvals,Lnnz,permout,nnzlplus,DIAG,LDL_D);

#Using  all input as output would not work because of unknown Lnnz before the factorization!!!
ccall((:selinv2julia,"/home/grad/wgong/Spatial/newSelInv3/SelInvCopy/EXAMPLES/julia.so"),
      void, (Int32,Int32,Ref{Cint},Ref{Cint},Ref{Cdouble},Ref{Cint},Ref{Cint},Ref{Cint},Ref{Cdouble},Ref{Cdouble}),
                        nnodes,nnz,colptr,rowind, nzvals,Lnnz,permout,nnzlplus,DIAG,LDL_D);
=#
VV= ccall((:selinv2julia,"/home/marcin/pyMRA/julia/SelInv/EXAMPLES/julia.so"),
        Ptr{Cdouble}, (Int32,Int32,Ref{Cint},Ref{Cint},Ref{Cdouble},Ref{Cint},Ref{Cint},Ref{Cint},Ref{Cdouble},Ref{Cdouble}),
                        nnodes,nnz,colptr,rowind, nzvals,Lnnz,permout,nnzlplus,DIAG,LDL_D);
v2 = pointer_to_array(VV, (3*Lnnz[1]+nnzlplus[1]))# for julia-0.4.5
## call C function to free memory allocated in C
#v2=unsafe_wrap(Array, VV, (3*Lnnz[1]+nnzlplus[1]),true);# was pointer_to_array, now deprecated in julia-0.5
#This unsafe_wrap function is labelled "unsafe" because it will crash if `pointer` is not a valid memory address to data of the requested length.
nnzInv = Lnnz[1];
iInv = v2[1:nnzInv];
jInv = v2[(nnzInv+1):(2*nnzInv)];
iInv =  convert(Array{Int32,1},iInv);
jInv = convert(Array{Int32,1},jInv);
LDL_L =  v2[(2*nnzInv+1):(3*nnzInv)];
invElement =  v2[(3*nnzInv+1):(3*nnzInv+nnzlplus[1])];
inv_pos = sparse(iInv,jInv,1);
k=Int64[1]
newLNZ = zeros(nnzInv);
for i=1:nnzInv
    if i == 1
       k = convert(Int64,DIAG[jInv[i]]);
    elseif jInv[i]>jInv[i-1]
       k = convert(Int64,DIAG[jInv[i]]);
    end
    newLNZ[i] = invElement[k];
    k = k + 1;
end

for i=1:nnzInv
    iInv[i] = permout[iInv[i]];
    jInv[i] = permout[jInv[i]];
end

## get the LDL component
#LDL_L_matrix = sparse(iInv,jInv,LDL_L);
#for i=1:nnodes
#    LDL_L_matrix[i,i] = 1;
#end
#LDL_D_matrix = sparse([1:1:nnodes],[1:1:nnodes],LDL_D);

##log determinent
logDet=sum(log(LDL_D))
#logDet =sum(log(diag((LDL_D_matrix))))

## reconstruct testM
#testM_reconstruct = LDL_L_matrix*LDL_D_matrix*transpose(LDL_L_matrix);#LDLT
#testM_reconstruct-testM_sparse

testM_inv = sparse(iInv,jInv,newLNZ);
testM_inv = testM_inv + transpose(testM_inv);
for i=1:nnodes
    testM_inv[i,i] = testM_inv[i,i]/2
end

## cholesky decomposition from LDLT
#cholLDL=Any
#cholLDL=(LDL_D_matrix.^0.5)*LDL_L_matrix'
ccall((:freeV,"/home/grad/wgong/Spatial/newSelInv3/SelInvCopy/EXAMPLES/julia.so"),
        Void, (Ptr{Cdouble},),VV) ## to free the space allocated in C

  return(testM_inv,
          logDet
         )

end


###################################################################################################3
using Distances  #for distance computation (pairwise())
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
     ceil(Int64,x)#iceil(x)
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


## simulate stationary data on 1-D grid (based on function DLLoglikelihood in R package ltsa)
function dllik(r,z)
  n = length(z)
  error = zeros(n)
  sigmasq = zeros(n)
  error[1] = z[1]
  sigmasq[1] = r[1]
  phi=zeros(1)
  phi[1] = r[2]/r[1]
  error[2] = z[2] - phi[1] * z[1]
  sigmasqkm1 = r[1] * (1 - phi[1]^2)
  sigmasq[2] = sigmasqkm1
  logg = log(r[1]) + log(sigmasqkm1)
  for k=2:(n-1)
    phikk = (r[k + 1] .- phi' * r[k:-1:2])./sigmasqkm1
    sigmasqk = sigmasqkm1 .* (1 - phikk.^2)
    phi -= phikk .* reverse(phi)
    append!(phi,phikk)
    sigmasqkm1 = sigmasqk
    logg += log(sigmasqk)
    error[k + 1] = (z[k + 1] - phi'*z[k:-1:1])[1]
    sigmasq[k + 1] = sigmasqk[1]
  end
  S = sum((error .* error)./sigmasq)
  LogL = S + logg
  return LogL[1]
end


# distance for either 1D or 2D locations.
function dist(loc1,loc2)
  return sqrt(sum((loc1-loc2).^2))
end

#distance for each pair of locations, output dense matrix
function distmat_d(locs1,locs2)
  m=size(locs1,1)
  n=size(locs2,1)
  distmat = zeros(Float64,m,n)
  for i=1:m
    for j=1:n
      distmat[i,j] = dist(locs1[i],locs2[j])
    end
  end
  return distmat
end

# determine nonzero(=1) entries within certain distance,output sparse matrix with binary entries
function nzmatrix(loc1,loc2,max)

  loc1=sort(loc1)
  loc2=sort(loc2)
  m=size(loc1,1)
  n=size(loc2,1)

  colptr = Int64[]
  rowind = Int64[]
  nzval = Float64[]

  ip=1
  bi = 1
  for j = 1:n #column
    append!(colptr,[ip])
    #find new bi, loop starts from previous bi
    while (dist(loc1[bi],loc2[j]) > max) & (bi < m)
    bi = bi + 1
    end
    for i = bi:m #row
      d = dist(loc1[i],loc2[j])
      #if d==0 d=e0 end
      if d <= max
        append!(rowind,[i])
        ip += 1
        append!(nzval,[1])
      else
        break
      end
    end
  end
  append!(colptr,[length(rowind)+1])
  return  SparseMatrixCSC(m, n, colptr,rowind , nzval)
end



#e0=1e-100
#e0/0.05<1e-40 #true
#=loc1=NewKnots[1]
loc2=NewKnots[1]
max=Newranges[2]
loc1=knots[1]
loc2=knots[1]
max=ranges[end]
M=distmatrix(loc1,loc2,max)
spy(M)
NewKnots[2]-knots[end]
M=distmatrix(NewKnots[2],NewKnots[2],Newranges[2])
M=distmatrix00(knots[1],knots[1],Newranges[2])
spy(M)
=#
# calculate distance within certain limit(max) for 1D locs, output sparse distance matrix.
function distmatrix(loc1,loc2,max)
  loc1=sort(loc1)
  loc2=sort(loc2)
  m=size(loc1,1)
  n=size(loc2,1)

  colptr = Int64[]
  rowind = Int64[]
  nzval = Float64[]

  ip=1
  bi = 1
  for j = 1:n #column
    append!(colptr,[ip])
    #find new bi, loop starts from previous bi
    while (euclidean(loc1[bi,:],loc2[j,:]) > max) & (bi < m)
    bi = bi + 1
    end
    for i = bi:m #row
      d = euclidean(loc1[i,:],loc2[j,:])
      if d==0 d=[1-1] end
      if d[1] <= max
        append!(rowind,[i])
        ip += 1
        append!(nzval,[d[1]])
      else
        break
      end
    end
  end
  append!(colptr,[length(rowind)+1])
  return  SparseMatrixCSC(m, n, colptr,rowind , nzval)
end
# calculate distance within certain limit(max) for 2D locs, output sparse distance matrix.
function distmatrix2D0(loc1,loc2,max)
  m=size(loc1,1)
  n=size(loc2,1)
  colptr = Int64[]
  rowind = Int64[]
  nzval = Float64[]
  ip=1
  bi = 1
  for j = 1:n #column
    append!(colptr,[ip])
    #find new bi, loop starts from previous bi
    #while ((euclidean(loc1[bi,1],loc2[j,1]) > max)  & (bi < m))
    #bi = bi + 1 # find the first row index that dist> max
    #end
    for i = 1:m #row
      #dx= euclidean(loc1[i,1],loc2[j,1])# dist on x axis
      #dy= euclidean(loc1[i,2],loc2[j,2])# dist on x axis
      d = euclidean(loc1[i,:],loc2[j,:])
      #if (dx>max)||(dy>max)
      #  break  #not work b/c break both loop
      #else
        if d==0 d=[1-1] end
        if d[1] <= max
          append!(rowind,[i])
          ip += 1
          append!(nzval,[d[1]])
        #else # if d > max doesn't mean that the following would also > max
        #  break
        end
      end
    end

  append!(colptr,[length(rowind)+1])
  return  SparseMatrixCSC(m, n, colptr,rowind , nzval)
  end


#function for tapering with max distance, for 2D domain
function distmatrix2D(loc1,loc2,max)
  loc1=sortrows(loc1,by=x->(x[1],x[2]))
  loc2=sortrows(loc2,by=x->(x[1],x[2]))
  m=size(loc1,1)
  n=size(loc2,1)
  colptr = Int64[]
  rowind = Int64[]
  nzval = Float64[]
  ip=1
  bi = 1
  for j = 1:n #column
    append!(colptr,[ip])
    #find new bi, loop starts from previous bi
    while ((euclidean(loc1[bi,1],loc2[j,1]) > max)  & (bi < m))
    bi = bi + 1 # find the first row index that dist> max
    end
    for i = bi:m #row
      #dx= euclidean(loc1[i,1],loc2[j,1])# dist on x axis
      #dy= euclidean(loc1[i,2],loc2[j,2])# dist on x axis
      d = euclidean(loc1[i,:],loc2[j,:])
      #if (dx>max)||(dy>max)
      #  break  #not work b/c break both loop
      #else
        if d==0 d=[1-1] end
        if d[1] <= max
          append!(rowind,[i])
          ip += 1
          append!(nzval,[d[1]])
        #else # if d > max doesn't mean that the following would also > max
        #  break
        end
      end
    end

  append!(colptr,[length(rowind)+1])
  return  SparseMatrixCSC(m, n, colptr,rowind , nzval)
end
#function for tapering with max distance, for 2D domain
function distmatrix2D01(loc1,loc2,max)
  loc1=sortrows(loc1,by=x->(x[1],x[2]))
  loc2=sortrows(loc2,by=x->(x[1],x[2]))
  m=size(loc1,1)
  n=size(loc2,1)
  colptr = Int64[]
  rowind = Int64[]
  nzval = Float64[]
  ip=1
  bi = 1
  for j = 1:n #column
    append!(colptr,[ip])
    #find new bi, loop starts from previous bi
    while ((euclidean(loc1[bi,1],loc2[j,1]) > max)  & (bi < m))
    bi = bi + 1 # find the first row index that dist> max
    end
    for i = bi:m #row
      #dx= euclidean(loc1[i,1],loc2[j,1])# dist on x axis
      #dy= euclidean(loc1[i,2],loc2[j,2])# dist on x axis
      d = euclidean(loc1[i,:],loc2[j,:])
      #if (dx>max)||(dy>max)
      #  break  #not work b/c break both loop
      #else
        if d==0 d=[1-1] end
        if d[1] <= max
          append!(rowind,[i])
          ip += 1
          append!(nzval,1)
        #else # if d > max doesn't mean that the following would also > max
        #  break
        end
      end
    end

  append!(colptr,[length(rowind)+1])
  return  SparseMatrixCSC(m, n, colptr,rowind , nzval)
end

## selected inversion of sparse spd matrix --- replace by fortran code https://math.berkeley.edu/~linlin/software.html
function selInv(A,B)
  inv(full(A)).*B  # B is matrix of ones and zeros that specifies the sparsity structure
end

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
#not work for 2D
function distfun(locs1,locs2)
  pairwise(Euclidean(),locs1',locs2')
end
#not work for pairwise
function distfun1(loc1,loc2)#was dist()
  return sqrt(sum((loc1-loc2).^2))
end
#
#locs1=[2,2]; locs2=[6,6]
#distfun(locs1,locs2)
#dist(locs1,locs2)

#covariance function for sparse distance matrix.
function co_spmat(dists,sig,prange,sm,sigeps)
  m=size(dists,1)
  n=size(dists,2)
  return  SparseMatrixCSC(m, n, dists.colptr, dists.rowval, covFunWE(dists.nzval,sig,prange,sm,sigeps))
end;


# Kanter function for sparse distance matrix
function kanter_spmat(dists) # tapering function of Kanter (1997)
  m=size(dists,1)
  n=size(dists,2)
  t=dists.nzval
  twopit=2*pi*t + 1.0e-300;  # add tiny number to ensure kanter(0)=1
  R=(1-t).*sin(twopit)./twopit+(1/pi)*(1-cos(twopit))./twopit;
  R[t.>1]=0; R[t.==0]=1;
  return  SparseMatrixCSC(m, n, dists.colptr, dists.rowval, R)
end

#taper function for sparse matrix with range
function taper_spmat(distances,range)
  return kanter_spmat(distances/range)
end




# tapering function of Kanter (1997)
function kanter(t)
# valid in R^3, minimizes the curvature at the origin, only has 2 derivatives at the
#    origin (not useful for matern with nu>2), for more details, see gneiting (2002), eq. 22
  twopit=2*pi*t + 1.0e-300;  # add tiny number to ensure kanter(0)=1
  R=(1-t).*sin(twopit)./twopit+(1/pi)*(1-cos(twopit))./twopit;
  R[t.>1]=0; R[t.==0]=1;
  return R
end

#taper function of Kanter with range
function taper(distances,range)
  return kanter(distances/range)
end

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
################################################################################################
# prior function for MRAtaper
function MRAprior_taper(dmat,dmat1,dmat_SelInv,knots,ranges,sig,prange,sm,sigeps)#nzmat,nzmat1
  M=length(ranges)-1
  Plist=Array(SparseMatrixCSC,M+1);  #a collection of sparse matrices
  LogDet=Array(Float64,M+1)
  K=Array(SparseMatrixCSC,M+1)
  V=deepcopy(dmat)#nzmat #deepcopy
  A=Any
  Vt=deepcopy(dmat1)#nzmat1 #deepcopy
  #Pchol=Array(Any,M+1)
  for m=0:M
     for l=0:m
        dists=dmat[m+1][l+1]
        V[m+1][l+1] =co_spmat(dists,sig,prange,sm,sigeps).*taper_spmat(dists,ranges[1])
        #dists=distfun(knots[m+1],knots[l+1])
        #V[m+1][l+1]=sparse(co(knots[m+1],knots[l+1]).*taper(dists,ranges[1]))
       for k=1:l
          A=(K[k]*Vt[m+1][k])######
          Vrows = (V[m+1][l+1]).rowval #updated whenever Vals changed
          Vvals = (V[m+1][l+1]).nzval #updated whenever Vals changed
          for j=1:size(V[m+1][l+1],2)
            i=nzrange(V[m+1][l+1],j) #return the range of indices of nz of a sp mat column
            dists1=distfun(knots[m+1][Vrows[i],:],knots[l+1][j,:]) #knots[l+1][j,:]' for j0.5
            Vvals[i] =   Vvals[i]- sum(broadcast(.*,A[:,Vrows[i]],Vt[l+1][k][:,j]),1)' #Vvals[i] -= (A[:,Vrows[i]]'*Vt[l+1][k][:,j])#################
            Vvals[i] =   Vvals[i] .* taper(dists1,ranges[k+1])####################Multiply the ratio from taper function, sparsity pattern already got from V
          end
        end
        Vt[m+1][l+1] = V[m+1][l+1]'
    end
    Plist[m+1]=.5*V[m+1][m+1]+.5*V[m+1][m+1]';
    (K[m+1],LogDet[m+1])=SIplusD(Plist[m+1],dmat_SelInv[m+1]);
  #  K[m+1]=sparse(inv(full(Plist[m+1])));
  #  LogDet[m+1]=logdet(full(Plist[m+1]))
  end

  #Stack Basis matrices and  precision matrices
  B=V[M+1][1]
  for m=1:M
    B=hcat(B,V[M+1][m+1])
  end
  P=Plist[1] #sparse(Pchol[1])*sparse(Pchol[1])'
  for m=1:M
    P=blkdiag(P,Plist[m+1]) #sparse(Pchol[m+1])*sparse(Pchol[m+1])')
  end

  return (B,V,Vt,P,K,LogDet)
end

#
function MRAprior_taper1(dmat,dmat1,dmat_SelInv,knots,locs,ranges,sig,prange,sm,sigeps)#nzmat,nzmat1
  M=length(ranges)-1
  Plist=Array(SparseMatrixCSC,M+1);  #a collection of sparse matrices
  LogDet=Array(Float64,M+1)
  K=Array(SparseMatrixCSC,M+1)
  V=deepcopy(dmat)#nzmat #deepcopy
  A=Any
  Vt=deepcopy(dmat1)#nzmat1 #deepcopy
  #Pchol=Array(Any,M+1)
  for m=0:M
     for l=0:m
        dists=dmat[m+1][l+1]
        V[m+1][l+1] =co_spmat(dists,sig,prange,sm,sigeps).*taper_spmat(dists,ranges[1])
        #dists=distfun(knots[m+1],knots[l+1])
        #V[m+1][l+1]=sparse(co(knots[m+1],knots[l+1]).*taper(dists,ranges[1]))
       for k=1:l
          A=(K[k]*Vt[m+1][k])######
          Vrows = (V[m+1][l+1]).rowval #updated whenever Vals changed
          Vvals = (V[m+1][l+1]).nzval #updated whenever Vals changed
          for j=1:size(V[m+1][l+1],2)
            i=nzrange(V[m+1][l+1],j) #return the range of indices of nz of a sp mat column
            dists1=distfun(locs[Vrows[i],:],knots[l+1][j,:]) #knots[l+1][j,:]' for j0.5
            Vvals[i] =   Vvals[i]- sum(broadcast(.*,A[:,Vrows[i]],Vt[l+1][k][:,j]),1)' #Vvals[i] -= (A[:,Vrows[i]]'*Vt[l+1][k][:,j])#################
            Vvals[i] =   Vvals[i] .* taper(dists1,ranges[k+1])####################Multiply the ratio from taper function, sparsity pattern already got from V
          end
        end
        Vt[m+1][l+1] = V[m+1][l+1]'
    end
    Plist[m+1]=.5*V[m+1][m+1]+.5*V[m+1][m+1]';
    (K[m+1],LogDet[m+1])=SIplusD(Plist[m+1],dmat_SelInv[m+1]);
  #  K[m+1]=sparse(inv(full(Plist[m+1])));
  #  LogDet[m+1]=logdet(full(Plist[m+1]))
  end

  #Stack Basis matrices and  precision matrices
  B=V[M+1][1]
  for m=1:M
    B=hcat(B,V[M+1][m+1])
  end
  P=Plist[1] #sparse(Pchol[1])*sparse(Pchol[1])'
  for m=1:M
    P=blkdiag(P,Plist[m+1]) #sparse(Pchol[m+1])*sparse(Pchol[m+1])')
  end

  return (B)
end

   #=
   function MRAprior0(co_spmat,taper_spmat,taper,dmat,dmat1,dmat_SelInv,knots,ranges)#nzmat,nzmat1
   Plist=Array(SparseMatrixCSC,M+1);  #a collection of sparse matrices
   Pchol=Array(Any,M+1)
   K=Array(SparseMatrixCSC,M+1)
   LogDet=Array(Float64,M+1)
   V=deepcopy(dmat)
   A=Any
   Vt=deepcopy(dmat1)
   for m=0:M
      for l=0:m
         dists=dmat[m+1][l+1]
         V[m+1][l+1]=co_spmat(dists).*taper_spmat(dists,ranges[1])
        for k=1:l
           A=(K[k]*Vt[m+1][k])
           Vrows = (V[m+1][l+1]).rowval #updated whenever Vals changed
           Vvals = (V[m+1][l+1]).nzval #updated whenever Vals changed
           for j=1:size(V[m+1][l+1],2)
             i=nzrange(V[m+1][l+1],j) #return the range of indices of nz of a sp mat column
            Vvals[i] -= (A[:,Vrows[i]]'*Vt[l+1][k][:,j])
             Vvals[i] .*= taper(dists[Vrows[i],j],ranges[k+1])
           end
         end
         Vt[m+1][l+1] = V[m+1][l+1]'
     end
     Plist[m+1]=.5*V[m+1][m+1]+.5*V[m+1][m+1]';
     #Pchol[m+1]=cholfact(Plist[m+1])#sparsechol(Plist[m+1])#
     #if m<M;
       (K[m+1],LogDet[m+1]) =SIplusD(Plist[m+1],dmat_SelInv[m+1]);
     #end
   end

   B=V[M+1][1]
   for m=1:M
     B=hcat(B,V[M+1][m+1])
   end
   P=Plist[1] #sparse(Pchol[1])*sparse(Pchol[1])'
   #spy(Plist[4])
   for m=1:M
     P=blkdiag(P,Plist[m+1]) #sparse(Pchol[m+1])*sparse(Pchol[m+1])')
     end

     return (B,V,P,K,LogDet)
   end
   =#
# calculate taper ranges given # of non-zeros per row
   function TaperRanges(nzrowprior,r0,J,M,n)
     #nzrowprior=2; # desired nz elements per row in prior
     rangemultiplier=nzrowprior/r0/2;
     range0=(domain[2]-domain[1])*rangemultiplier;
     #=
     if M==1
       #rs=[r0*M; n;];
       ranges=[range0; range0/J^M];
     else
     =#
       #rs=r0*J.^[0:M;]; #rs was used for generate knots, not needed now
       ranges=range0./J.^[0:M;]
     #end
     #ranges=1e100*ones(M+1) #to test no tapering
     return ranges
   end
# calculate taper ranges based on effective range
   function TaperRangesNew(sig,prange,sm,r0,J,M,k1)
     ERd=0.05
     ER=effrange(sig,prange,sm)
     ## find smallest k s.t. ranges0=k*ER
     ##k=4
     range0=ER*k1
     ranges=range0./J.^[0:M;]
     return ranges
   end
# calculate taper ranges with given effective range
   function TaperRangesNew0(ER,r0,J,M,k1)
     #ERd=0.05
     #ER=effrange(sig,prange,sm)
     ## find smallest k s.t. ranges0=k*ER
     ##k=4
     range0=ER*k1
     ranges=range0./J.^[0:M;]
     return ranges
   end
#=
   ##ranges compute according to # of knots each res.
   TaperRangesNew1=function(NewKnots,nzrowprior)
     ranges=zeros(2)
     rangemultiplier=nzrowprior/2;
     ranges[1]=(domain[2]-domain[1])*rangemultiplier/length(NewKnots[1]);
     ranges[2]=(domain[2]-domain[1])*rangemultiplier/length(NewKnots[2])
   return ranges
   end
=#
# sparsity for MRA taper prior
   function TaperMat(distmatrix,locs,knots,ranges,J,M)
     spbdS=Array(Any,M+1)
     spbdS1=Array(Any,M+1)
    # if sigeps==0
       #spbS=Array(Any,M+1)
       #spbdS=Array(Any,M+1)
       #spbdS1=Array(Any,M+1)

       for m=0:M
         #spbS[m+1]=distmatrix(locs,knots[m+1],ranges[m+1])
         spbdS[m+1]=distmatrix(locs,knots[m+1],ranges[m+1])
         spbdS1[m+1]=distmatrix(knots[m+1],locs,ranges[m+1])
       end

    # end
     #################   calculate prior for nz elements    ###################
     # pre-calculate sparse matrices for SelInv
     dmat_SelInv=Array(Any,M+1)
     for m=0:M
       dmat_SelInv[m+1]=distmatrix(knots[m+1],knots[m+1],(2+2/J)*ranges[m+1])
     end

     # pre-calculate sparse distance matrices
     dmat=Array(Any,M+1)
     for m=0:M
       dmat[m+1]=Array(SparseMatrixCSC,m+1)
       for l=0:m #########################loop for knots btw resolutions.
         dmat[m+1][l+1]=distmatrix(knots[m+1],knots[l+1],ranges[l+1])
       end
     end
     ##########Switch row and col knots
     dmat1=Array(Any,M+1)
     for m=0:M
       dmat1[m+1]=Array(SparseMatrixCSC,m+1)
       for l=0:m #########################loop for knots btw resolutions.
         dmat1[m+1][l+1]=distmatrix(knots[l+1],knots[m+1],ranges[l+1])
       end
     end

     return (dmat,dmat1,dmat_SelInv,spbdS,spbdS1)
   end
#=
#mannully make the dist matrix symmetric
   function TaperMat00(locs,knots,ranges,J,M,sigeps)
     spbdS=Array(Any,M+1)
     spbdS1=Array(Any,M+1)
     if sigeps==0
       #spbS=Array(Any,M+1)
       #spbdS=Array(Any,M+1)
       #spbdS1=Array(Any,M+1)
       for m=0:M-1
         #spbS[m+1]=distmatrix(locs,knots[m+1],ranges[m+1])
         spbdS[m+1]=distmatrix(locs,knots[m+1],ranges[m+1])
         spbdS1[m+1]=distmatrix(knots[m+1],locs,ranges[m+1])
       end
       spbdS[M+1]=distmatrix00(locs,knots[M+1],ranges[M+1])
       spbdS1[M+1]=distmatrix00(knots[M+1],locs,ranges[M+1])

     end
     #################   calculate prior for nz elements    ###################
     # pre-calculate sparse matrices for SelInv
     dmat_SelInv=Array(Any,M+1)
     for m=0:M
       dmat_SelInv[m+1]=distmatrix00(knots[m+1],knots[m+1],(2+2/J)*ranges[m+1])
     end

     # pre-calculate sparse distance matrices
     dmat=Array(Any,M+1)
     for m=0:M
       dmat[m+1]=Array(SparseMatrixCSC,m+1)
       for l=0:m-1 #########################loop for knots btw resolutions.
         dmat[m+1][l+1]=distmatrix(knots[m+1],knots[l+1],ranges[l+1])
       end
       dmat[m+1][m+1]=distmatrix00(knots[m+1],knots[m+1],ranges[m+1])
     end
     ##########Switch row and col knots
     dmat1=Array(Any,M+1)
     for m=0:M
       dmat1[m+1]=Array(SparseMatrixCSC,m+1)
       for l=0:m-1 #########################loop for knots btw resolutions.
         dmat1[m+1][l+1]=distmatrix(knots[l+1],knots[m+1],ranges[l+1])
       end
        dmat1[m+1][m+1]=distmatrix00(knots[m+1],knots[m+1],ranges[m+1])
     end
     return (dmat,dmat1,dmat_SelInv,spbdS,spbdS1)
   end

=#
## logliklihood of MRA taper with measurement error for 1D and 2D domain
function MRAtaperWE(ranges,knots,locs,z,r0,J,M,sig,prange,sm,sigeps)#distmatrix,MRAprior_taper,MRApriorS,MRAloglik,
   n=length(z)
  if size(locs,2)==1
     (dmat,dmat1,dmat_SelInv)=TaperMat(distmatrix,locs,knots,ranges,J,M,sigeps)[1:3]
     TNNZ=TotNNZ1(dmat_SelInv)
    #compute prior
     tic()
     (B,V,Vt,P,K,LogDet)=MRAprior_taper(dmat,dmat1,dmat_SelInv,knots,ranges,sig,prange,sm,sigeps)#nzmat,nzmat1 # output K for prediction
     T=toq()
     #compute logliklihood of MRAtaper
     loglik_MRTA=MRAloglikWE(P,B,LogDet,sigeps,z)
     loglikMRTA=0
elseif size(locs,2)==2
     (dmat,dmat1,dmat_SelInv)=TaperMat(distmatrix2D,locs,knots,ranges,J,M,sigeps)[1:3]
     TNNZ=TotNNZ1(dmat_SelInv)
    #compute prior
     tic()
     (B,V,Vt,P,K,LogDet)=MRAprior_taper(dmat,dmat1,dmat_SelInv,knots,ranges,sig,prange,sm,sigeps)#nzmat,nzmat1 # output K for prediction
     T=toq()
     #compute logliklihood of MRAtaper
     loglik_MRTA=MRAloglikWE(P,B,LogDet,sigeps,z)
     loglikMRTA=0
  end
    return (loglikMRTA,loglik_MRTA,T,TNNZ)#,ranges[end])
end
