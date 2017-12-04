########################################################################
### MRAblock/1RAblock estimation, prediction on SST data application
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
function estMB(ir0,set) ## ir0 to iterate over different settings of MRA; set denotes datasets for repetation
  mod="MRAblock" ## mod: "1RAblock" or "MRAblock"
## initiate parameters
  sm=0.5 # smoothness
  theta=[0.0466426,0.340289,4.36922e-11]
  sig=theta[1] # var
  prange=theta[2]# range parameter
  sigeps=theta[3] # nugget
## different settings of MRA for comparison
  r0s=[3,6,5,10,11,12]
  Ms=[5,5,5,4,4,4]
  Js=[2,2,2,2,2,2]
## settings for 1RA
  #r0s=[5,10,10,10,20]
  #Ms=[1,1,1,1,1]
  #Js=[30,15,8,6,15] ## approximately sqrt(n/r0^2)
  iJ=iM=ir0
  r0=r0s[ir0]
  M=Ms[iM]
  J=Js[iJ]

## read data and generate knots
  datalocs=readcsv("/home/grad/wgong/Spatial/Data/SSTmodel"*string(set)*".csv")[2:end,2:4]#35456
  z=collect(Float64,datalocs[:,3])
  n=length(z)
  locs=convert(Array{Float64,2}, datalocs[:,1:2])
  domainX=[minimum(locs[:,1]);maximum(locs[:,1]);];
  domainY=[minimum(locs[:,2]);maximum(locs[:,2]);];
  knots=D2knots(r0,J,M,domainX,domainY)
  knots[M+1]=deepcopy(locs)
  (boundsX,boundsY)= Bounds2D(domainX,domainY,J,M)
## sparsity structure for MRA
  (spbd,spbd1,sp_Selinv,spbdS,spbdS1)= BlockMat(spBd,knots,locs,boundsX,boundsY)
## MRA loglikelihood with initial parameters for compilation purpose
  tic()
  llb=MRAblockll(knots,locs,z,r0,J,M,theta,sm,spbd,spbd1,sp_Selinv,spbdS,spbdS1)
  timeB=toq()
################################ parameters estimation  ########################
## inital parameters
  sm=0.5
  theta=[0.0466426,0.340289,4.36922e-11]
  init=log(theta)
## optimization
  tic();
  optimresB=optimize(x-> nloglikB(knots,locs,z,r0,J,M,x,sm,spbd,spbd1,sp_Selinv,spbdS,spbdS1),
                      g_tol = 1e-5,init,iterations=500,store_trace=true,show_trace=true,show_every=1)
  totaltimeB=toc();
  loglikelihoodB= -optimresB.minimum
  theta_hatB =exp(optimresB.minimizer)
## write setting and estimations to file
  filePath = "/home/grad/wgong/Spatial/Data/SST_estMB"*string(set)*"_myid"
  #filePath = "/home/grad/wgong/Spatial/Data/SST_est1B"*string(set)*"_myid" ## 1B for 1RAblock
  fid1 =  open(filePath*string(myid())*".csv", "a")
  writecsv( fid1, hcat(r0,J,M,n,sm,llb,timeB,loglikelihoodB,totaltimeB,
                        init[1],init[2],init[3],theta_hatB[1],theta_hatB[2],theta_hatB[3],mod,set)  )
  close(fid1)
end #function estMB

#######################################################
## MRA prediction ##
#######################################################
function predMB(ir0,set)
  mod="MRAblock" ## mod: "1RAblock" or "MRAblock"
## initiate parameters
  sm=0.5 # smoothness
  theta=[0.0466426,0.340289,4.36922e-11]
  sig=theta[1] # var
  prange=theta[2]# range parameter
  sigeps=theta[3] # nugget
## different settings of MRA for comparison
  r0s=[3,6,5,10,11,12]
  Ms=[5,5,5,4,4,4]
  Js=[2,2,2,2,2,2]
## settings for 1RA
  #r0s=[5,10,10,10,20]
  #Ms=[1,1,1,1,1]
  #Js=[30,15,8,6,15] ## approximately sqrt(n/r0^2)
  iJ=iM=ir0
  r0=r0s[ir0]
  M=Ms[iM]
  J=Js[iJ]
## read data and generate knots
  datalocs=readcsv("/home/grad/wgong/Spatial/Data/SSTmodel"*string(set)*".csv")[2:end,2:4]#35456
  z=collect(Float64,datalocs[:,3])
  n=length(z)
  locs=convert(Array{Float64,2}, datalocs[:,1:2])
  domainX=[minimum(locs[:,1]);maximum(locs[:,1]);];
  domainY=[minimum(locs[:,2]);maximum(locs[:,2]);];
  knots=D2knots(r0,J,M,domainX,domainY)
  knots[M+1]=deepcopy(locs)
  (boundsX,boundsY)= Bounds2D(domainX,domainY,J,M)
## sparsity structure for MRA
  (spbd,spbd1,sp_Selinv,spbdS,spbdS1)= BlockMat(spBd,knots,locs,boundsX,boundsY)
## read theta_hat from the file contains estimtes
  est=DataFrame(readcsv("/home/grad/wgong/Spatial/Data/catSST_estMB"*string(set)*"all.csv"))# or est1B for 1RA
  #est=DataFrame(readcsv("/home/grad/wgong/Spatial/Data/catSST_est1B"*string(set)*"all.csv"))
  rename!(est,names(est),[:r0,:J,:M,:n,:sm,:llb,:timeB,:loglikelihoodB,:totaltimeB,
                     :init1,:init2,:init3,:theta_hatB1,:theta_hatB2,:theta_hatB3,:MRAtype,:set])
  theta_hatB=Array(est[(est[:r0].==r0)&(est[:J].==J)&(est[:M].==M)&(est[:MRAtype].==mod),
                    [:theta_hatB1,:theta_hatB2,:theta_hatB3,:loglikelihoodB]][1,:])'

####################################### Prediction #################################
### Prediction with areal test data ###
 preddatalocs=readcsv("/home/grad/wgong/Spatial/Data/SSTpred_area"*string(set)*".csv")[2:end,2:4]#46579
 predlocs=convert(Array{Float64,2}, preddatalocs[:,1:2])
 preddata=collect(Float64,preddatalocs[:,3])
## sparsity structure for prediction
 (spbdS_p,spbdS1_p)=BlockMat_p(spBd,knots,predlocs,boundsX,boundsY)
## prediction with try and catch to adjust estimate slightly in case of numerical problem
 predresult= try pred_block(z,predlocs,theta_hatB,knots,sm,boundsX,boundsY,spbd,spbd1,sp_Selinv,spbdS,spbdS1,spbdS_p,spbdS1_p)
 catch err
    println(err)
    fac=log10(1e-12/theta_hatB[3])
    theta_hatB[3]=theta_hatB[3]*10^fac ## to increase the nugget estimate to 1e-12 for numerical stability
    predresult= pred_block(z,predlocs,theta_hatB,knots,sm,boundsX,boundsY,spbd,spbd1,sp_Selinv,spbdS,spbdS1,spbdS_p,spbdS1_p)
  end
## prediction
 tic()
 predresult=pred_block(z,predlocs,theta_hatB,knots,sm,boundsX,boundsY,spbd,spbd1,sp_Selinv,spbdS,spbdS1,spbdS_p,spbdS1_p)
 predtime=toq()
## write prediction results to file
  writecsv("/home/grad/wgong/Spatial/Data/newSSTpredresult_block_area"*string(r0,M)*"_set"*string(set)*".csv",predresult)
 #writecsv("/home/grad/wgong/Spatial/Data/newSSTpredresult_1block_area"*string(r0,J,M)*"_set"*string(set)*".csv",predresult) # 1block for 1RA

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
 filePath = "/home/grad/wgong/Spatial/Data/newSSTpredsum_block_area"*string(r0,M)*"_set"*string(set)*".csv" # 1block for 1RA
 #filePath = "/home/grad/wgong/Spatial/Data/newSSTpredsum_1block_area"*string(r0,J,M)*"_set"*string(set)*".csv" # 1block for 1RA

 fid1 = open(filePath, "a")
 writecsv(fid1,hcat(r0,J,M,n,sm,0,theta_hatB[4], #use range0=0 for 1block to match taper ranges
                    predtime,MSPE,RSPE,logscore,crps,mod,"area",set)) ## mod: "1RAblock" or "MRAblock"
 close(fid1)

#########################################################
## prediction for random location test data
 preddatalocs=readcsv("/home/grad/wgong/Spatial/Data/SSTpred_random"*string(set)*".csv")[2:end,2:4]
 predlocs=convert(Array{Float64,2}, preddatalocs[:,1:2])
 preddata=collect(Float64,preddatalocs[:,3])
## sparsity structure for prediction
 (spbdS_p,spbdS1_p)=BlockMat_p(spBd,knots,predlocs,boundsX,boundsY)
## prediction with try and catch to adjust estimate slightly in case of numerical problem
 predresult= try pred_block(z,predlocs,theta_hatB,knots,sm,boundsX,boundsY,spbd,spbd1,sp_Selinv,spbdS,spbdS1,spbdS_p,spbdS1_p)
 catch err
    println(err)
    fac=log10(1e-12/theta_hatB[3])
    theta_hatB[3]=theta_hatB[3]*10^fac ## to increase the nugget estimate to 1e-12 for numerical stability
    predresult= pred_block(z,predlocs,theta_hatB,knots,sm,boundsX,boundsY,spbd,spbd1,sp_Selinv,spbdS,spbdS1,spbdS_p,spbdS1_p)
  end
## prediction
 tic()
 predresult=pred_block(z,predlocs,theta_hatB,knots,sm,boundsX,boundsY,spbd,spbd1,sp_Selinv,spbdS,spbdS1,spbdS_p,spbdS1_p)
 predtime=toq()
## write prediction results to file
 writecsv("/home/grad/wgong/Spatial/Data/newSSTpredresult_block_random"*string(r0,M)*"_set"*string(set)*".csv",predresult)
 #writecsv("/home/grad/wgong/Spatial/Data/newSSTpredresult_1block_random"*string(r0,J,M)*"_set"*string(set)*".csv",predresult)

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
 filePath = "/home/grad/wgong/Spatial/Data/newSSTpredsum_block_random"*string(r0,M)*"_set"*string(set)*".csv" # 1block for 1RA
 #filePath = "/home/grad/wgong/Spatial/Data/newSSTpredsum_1block_random"*string(r0,J,M)*"_set"*string(set)*".csv"

 fid1 = open(filePath, "a")
  writecsv(fid1,hcat(r0,J,M,n,sm,0,theta_hatB[4],predtime,MSPE,RSPE,logscore, crps, ##range0=0 for 1block to match taper ranges
                     mod,"random",set)) ## mod: "1RAblock" or "MRAblock"
  close(fid1)
end #function predMB
