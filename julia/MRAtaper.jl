########################################################################
### MRAtaper/1RAtaper estimation, prediction on SST data application
########################################################################
##load packages and functions
using Optim, DataFrames, Distributions
include("/home/grad/wgong/Spatial/MRAcode/MRAfunctions.jl")
#include("/home/grad/wgong/Spatial/MRAcode/MRAfcts_taper_sig.jl")
#include("/home/grad/wgong/Spatial/MRAcode/sim_functions_sig.jl")
#include("/home/grad/wgong/Spatial/MRAcode/BlockFunctions_sig.jl")
#include("/home/grad/wgong/Spatial/MRAcode/est_fct_sig.jl")
############################################
## parameters estimation  ##
#############################################
function estMT(ir0,set) ## ir0 to iterate over different settings of MRA; set denotes datasets for repetation
  mod="MRAtaper" ## mod: "1RAtaper" or "MRAtaper"
## initiate parameters
  sm=0.5 # smoothness
  theta=[0.0466426,0.340289,4.36922e-11]
  sig=theta[1] # var
  prange=theta[2]# range parameter
  sigeps=theta[3] # nugget
## different settings of MRA for comparison
  r0s=[12,12,12,24,24,24]
  Ms=[4,4,4,3,3,3]
  Js=[2,2,2,2,2,2]
  ks=[2,3,4,1,2,3]
## settings for 1RA
  #r0s=[24,24,30,30,24,30]
  #Ms=[1,1,1,1,1,1,1]
  #Js=[8,8,8,15,8,8]
  #ks=[1,2,2,2,3,3]

  ik=iJ=iM=ir0
  r0=r0s[ir0]
  M=Ms[iM]
  J=Js[iJ]
  k=ks[ik]
  J2=J^2
## read data and generate knots
  datalocs=readcsv("/home/grad/wgong/Spatial/Data/SSTmodel"*string(set)*".csv")[2:end,2:4]#35456
  z=collect(Float64,datalocs[:,3])
  n=length(z)
  locs=convert(Array{Float64,2}, datalocs[:,1:2])
  domainX=[minimum(locs[:,1]);maximum(locs[:,1]);];
  domainY=[minimum(locs[:,2]);maximum(locs[:,2]);];
  knots=D2knots(r0,J,M,domainX,domainY)
  knots[M+1]=deepcopy(locs)

## sparsity structure for MRA
  range0=k
  ranges=range0./J2.^[0:M;]
  (dmat,dmat1,dmat_SelInv,spbdS,spbdS1)=TaperMat(distmatrix2D0,locs,knots,ranges,J,M)
## MRA loglikelihood with initial parameters for compilation purpose
  tic()
  llt=MRAtaperll(knots,locs,z,r0,J,M,ranges,theta,sm,dmat,dmat1,dmat_SelInv,spbdS,spbdS1)
  timeT=toq()
################################ parameters estimation  ########################
## inital parameters
  sm=0.5
  theta=[0.0466426,0.340289,4.36922e-11]
  init=log(theta)
## optimization
  tic();
  optimresT=optimize(x -> nloglikT(knots,locs,z,r0,J,M,ranges,x,sm,dmat,dmat1,dmat_SelInv,spbdS,spbdS1),
                      g_tol = 1e-5,init,iterations=500,show_trace=true,show_every=1)#store_trace=true,show_trace=true
  totaltimeT=toc();
  loglikelihoodT= -optimresT.minimum
  theta_hatT =exp(optimresT.minimizer)
## write setting and estimations to file
  filePath = "/home/grad/wgong/Spatial/Data/SST_estMT"*string(set)*"_myid"
  #filePath = "/home/grad/wgong/Spatial/Data/SST_est1T"*string(set)*"_myid" ## 1T for 1RAtaper
  fid1 =  open(filePath*string(myid())*".csv", "a")
  writecsv( fid1,hcat(r0,J,M,n,sm,llt,timeT,range0,loglikelihoodT,totaltimeT,
                      init[1],init[2],init[3],theta_hatT[1],theta_hatT[2],theta_hatT[3],mod,set) )
  close(fid1)
end #function estMT

#######################################################
## MRA prediction ##
#######################################################
function predMT(ir0,set)
  mod="MRAtaper" ## mod: "1RAtaper" or "MRAtaper"
## initiate parameters
  sm=0.5 # smoothness
  theta=[0.0466426,0.340289,4.36922e-11]
  sig=theta[1] # var
  prange=theta[2]# range parameter
  sigeps=theta[3] # nugget
## different settings of MRA for comparison
  r0s=[12,12,12,24,24,24]
  Ms=[4,4,4,3,3,3]
  Js=[2,2,2,2,2,2]
  ks=[2,3,4,1,2,3]
## settings for 1RA
  #r0s=[24,24,30,30,24,30]
  #Js=[8,8,8,15,8,8]
  #Ms=[1,1,1,1,1,1,1]
  #ks=[1,2,2,2,3,3]
  iJ=iM=ir0
  r0=r0s[ir0]
  M=Ms[iM]
  J=Js[iJ]
  J2=J^2
## read data and generate knots
  datalocs=readcsv("/home/grad/wgong/Spatial/Data/SSTmodel"*string(set)*".csv")[2:end,2:4]#35456
  z=collect(Float64,datalocs[:,3])
  n=length(z)
  locs=convert(Array{Float64,2}, datalocs[:,1:2])
  domainX=[minimum(locs[:,1]);maximum(locs[:,1]);];
  domainY=[minimum(locs[:,2]);maximum(locs[:,2]);];
  knots=D2knots(r0,J,M,domainX,domainY)
  knots[M+1]=deepcopy(locs)
## sparsity structure for MRA
  range0=k
  ranges=range0./J2.^[0:M;]
  (dmat,dmat1,dmat_SelInv,spbdS,spbdS1)=TaperMat(distmatrix2D0,locs,knots,ranges,J,M)
## read theta_hat from the file contains estimtes
  est=DataFrame(readcsv("/home/grad/wgong/Spatial/Data/catSST_estMT"*string(set)*"all.csv"))# or est1T for 1RA
  #est=DataFrame(readcsv("/home/grad/wgong/Spatial/Data/catSST_est1T"*string(set)*".csv"))
  rename!(est,names(est),[:r0,:J,:M,:n,:sm,:llt,:timeT,:range0,:loglikelihoodT,:totaltimeT,
                          :init1,:init2,:init3,:theta_hatT1,:theta_hatT2,:theta_hatT3,:MRAtype,:set])
  theta_hatT=Array(est[(est[:r0].==r0)&(est[:J].==J)&(est[:M].==M)&(est[:range0].==k)&(est[:MRAtype].==mod),
                    [:theta_hatT1,:theta_hatT2,:theta_hatT3,:loglikelihoodT]][1,:])'

####################################### Prediction #################################
### Prediction with areal test data ###
 preddatalocs=readcsv("/home/grad/wgong/Spatial/Data/SSTpred_area"*string(set)*".csv")[2:end,2:4]#46579
 predlocs=convert(Array{Float64,2}, preddatalocs[:,1:2])
 preddata=collect(Float64,preddatalocs[:,3])
## sparsity structure for prediction
 (spbdS_p,spbdS1_p)= TaperMat_p(distmatrix2D0,predlocs,knots,ranges)
## prediction with try and catch to adjust estimate slightly in case of numerical problem
  predresult= try pred_taper(distmatrix2D0,z,predlocs,theta_hatT,knots,ranges,sm,dmat,dmat1,dmat_SelInv,spbdS,spbdS1,spbdS_p,spbdS1_p)
  catch err
     println(err)
     fac=log10(1e-12/theta_hatT[3])
     theta_hatT[3]=theta_hatT[3]*10^fac ## to increase the nugget estimate to 1e-12 for numerical stability
     predresult= pred_taper(distmatrix2D0,z,predlocs,theta_hatT,knots,ranges,sm,dmat,dmat1,dmat_SelInv,spbdS,spbdS1,spbdS_p,spbdS1_p)
   end
## prediction
  tic()
  predresult=pred_taper(distmatrix2D0,z,predlocs,theta_hatT,knots,ranges,sm,dmat,dmat1,dmat_SelInv,spbdS,spbdS1,spbdS_p,spbdS1_p)
  predtime=toq()
## write prediction results to file
  writecsv("/home/grad/wgong/Spatial/Data/newSSTpredresult_taper_area"*string(r0,M)*"_set"*string(set)*".csv",predresult)
 #writecsv("/home/grad/wgong/Spatial/Data/try_newSSTpredresult_1taper_area"*string(r0,J,M,k)*"_set"*string(set)*".csv",predresult)

## calculate MSPE and RSPE
  MSPE=mean(preddata-predresult[:,3])^2
  RSPE=sqrt(mean(preddata-predresult[:,3])^2)
## calculate log-score, CRPS(continuous ranked probability score)
  predresult=DataFrame(hcat(predresult,preddata))
  rename!(predresult,  names(predresult)[3:5], [:mean,:sd,:z])
  llcusum=0
  crps_cusum=0
  for i=1:np
    sigma=predresult[i,:sd]
    mu=predresult[i,:mean]
    y=predresult[i,:z]
    mod=Normal(mu,sigma)
    crps_cusum = crps_cusum + sigma*( (y-mu)/sigma * (2*cdf(mod,y)-1)  + 2* pdf(mod,y)- 1/sqrt(pi) )
    llcusum = llcusum + logpdf(mod,y)
  end
  logscore=llcusum/np
  crps=crps_cusum/np
## write prediction summary to file
 filePath = "/home/grad/wgong/Spatial/Data/newSSTpredsum_taper_area"*string(r0,M)*"_set"*string(set)*".csv" # 1taper for 1RA
 #filePath = "/home/grad/wgong/Spatial/Data/try_newSSTpredsum_1taper_area"*string(r0,J,M,k)*"_set"*string(set)*".csv"

 fid1 = open(filePath, "a")
 writecsv(fid1,hcat(r0,J,M,n,sm,range0,theta_hatT[4],
                    predtime,MSPE,RSPE,logscore,crps,mod,"area",set)) ## mod: "1RAtaper" or "MRAtaper"
 close(fid1)

#########################################################
## prediction for random location test data
 preddatalocs=readcsv("/home/grad/wgong/Spatial/Data/SSTpred_random"*string(set)*".csv")[2:end,2:4]
 predlocs=convert(Array{Float64,2}, preddatalocs[:,1:2])
 preddata=collect(Float64,preddatalocs[:,3])
 ## sparsity structure for prediction
  (spbdS_p,spbdS1_p)= TaperMat_p(distmatrix2D0,predlocs,knots,ranges)
 ## prediction with try and catch to adjust estimate slightly in case of numerical problem
   predresult= try pred_taper(distmatrix2D0,z,predlocs,theta_hatT,knots,ranges,sm,dmat,dmat1,dmat_SelInv,spbdS,spbdS1,spbdS_p,spbdS1_p)
   catch err
      println(err)
      fac=log10(1e-12/theta_hatT[3])
      theta_hatT[3]=theta_hatT[3]*10^fac ## to increase the nugget estimate to 1e-12 for numerical stability
      predresult= pred_taper(distmatrix2D0,z,predlocs,theta_hatT,knots,ranges,sm,dmat,dmat1,dmat_SelInv,spbdS,spbdS1,spbdS_p,spbdS1_p)
    end
 ## prediction
   tic()
   predresult=pred_taper(distmatrix2D0,z,predlocs,theta_hatT,knots,ranges,sm,dmat,dmat1,dmat_SelInv,spbdS,spbdS1,spbdS_p,spbdS1_p)
   predtime=toq()
 ## write prediction results to file
 writecsv("/home/grad/wgong/Spatial/Data/newSSTpredresult_taper_random"*string(r0,M)*"_set"*string(set)*".csv",predresult)
# writecsv("/home/grad/wgong/Spatial/Data/try_newSSTpredresult_1taper_random"*string(r0,J,M,k)*"_set"*string(set)*".csv",predresult)

## calculate MSPE and RSPE
 MSPE=mean(preddata-predresult[:,3])^2
 RSPE=sqrt(mean(preddata-predresult[:,3])^2)
## calculate log-score, CRPS(continuous ranked probability score)
 predresult=DataFrame(hcat(predresult,preddata))
 rename!(predresult,  names(predresult)[3:5], [:mean,:sd,:z])
 llcusum=0
 crps_cusum=0
 for i=1:np
   sigma=predresult[i,:sd]
   mu=predresult[i,:mean]
   y=predresult[i,:z]
   mod=Normal(mu,sigma)
   crps_cusum = crps_cusum + sigma*( (y-mu)/sigma * (2*cdf(mod,y)-1)  + 2* pdf(mod,y)- 1/sqrt(pi) )
   llcusum = llcusum + logpdf(mod,y)
 end
 logscore=llcusum/np
 crps=crps_cusum/np
## write prediction summary to file
 filePath = "/home/grad/wgong/Spatial/Data/newSSTpredsum_taper_random"*string(r0,M)*"_set"*string(set)*".csv"
 #filePath = "/home/grad/wgong/Spatial/Data/try_newSSTpredsum_1taper_random"*string(r0,J,M,k)*"_set"*string(set)*".csv"

 fid1 = open(filePath, "a")
  writecsv(fid1,hcat(r0,J,M,n,sm,range0,theta_hatT[4],predtime,MSPE,RSPE,logscore, crps,
                     mod,"random",set)) ## mod: "1RAtaper" or "MRAtaper"
  close(fid1)
end #function predMT
