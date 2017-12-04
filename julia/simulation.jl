## loops for 2D simulation
function sim2D(ir0,seed)
r0s=[5,10,20]#160
Ms=[5,4,3]

#r0s=[3,6,12,4,8,16,5,10,20,3,6,12,24,7,14,28,2,4,8,16,32,9,18,36] #for 24, 32, 40, 48, 56, 64, 72
#Ms=[3,2,1,3,2,1,3,2,1,4,3,2,1,3,2,1,5,4,3,2,1,3,2,1]

#r0s=[7,14,28,8,16,32,9,18,36,10,20,40,11,22,44,6,12,24] # 112,128,144,160,176
#Ms=[4,3,2,4,3,2,4,3,2,4,3,2,4,3,2,4,3,2]

#r0s=[9,18,36,5,10,20,40,6,12,24,48]# 72 80 96
#Ms=[3,2,1,4,3,2,1,4,3,2,1]

#r0s=[3,6,12]#24
#Ms=[3,2,1]#

#r0s=[5,10,20,8,16,32,10,20,40,12,24,48]#80,128,160,192
#Ms=[4,3,2,4,3,2,4,3,2,4,3,2]

Js=2;
ks=2;
sms=[0.5,1.5,2.5];
pranges=[0.05,0.15];
    sigepss=[0.05,0.2]
filePath = "/home/marcin/pyMRA/julia/data/Exp_Theta0.1_X10_Y10.csv"
#filePath = "/home/grad/wgong/Spatial/TuningMRTA/parallel2/newD2_160nug02_$seed._"#112sm05_pr15$seed._"

for ism=1:length(sms),iprange=1:length(pranges),isigeps=1:length(sigepss),iJ=1:length(Js),ik=1:length(ks)
    iM=ir0
    sm=sms[ism]
    prange=pranges[iprange]
    sigeps=sigepss[isigeps]
    r0=r0s[ir0]
    M=Ms[iM]
    J=Js[iJ]
    k=ks[ik]
    J2=J^2
    n=r0*J^M;
    inputData = "/home/marcin/pyMRA/julia/data/Exp_Theta0.1_X10_Y10.csv"
    datalocs=readdlm(inputData,',')[2:end]
    z=collect(Float64,datalocs[:,4])
    print(z)
    locs=convert(Array{Float64,2}, datalocs[:,3])
    domain=[minimum(locs);maximum(locs);];
    loglikGP0=datalocs[:,4][1]
    if n<30000
    tic()
      GPll(co,locs,z,prange,sm,sigeps)
    TGP0=toq()
    else
      TGP0=0
    end
    knots=D2knots(r0,J,M,domain)
  # p=plotknots2D(knots[1:3])
###########################################################################
#MRA taper/block
  ranges=TaperRangesNew(prange,sm,r0,J2,M,k)
  (loglikMRTA,loglik_MRTA,TT,TNNZ) =  MRAtaperWE(ranges,knots,locs,z,r0,J2,M,prange,sm,sigeps)

  bounds= Bounds(domain,J,M)
  (loglik_MRBA,TB,BNNZ) = MRAblockWE(bounds,knots,locs,z,r0,J2,M,prange,sm,sigeps)
#1RA taper/block
###############################################################################
  (loglik_1RTA,T1RTA,TNNZT1 )= try OneRAtaper(oneRAknots0D2,domain,locs,z,prange, sm, sigeps, knots, r0, r0, ranges, 1,k)
  catch e
    if (isa(e,DomainError))
    (loglik_1RTA,T1RTA,TNNZT1 )= OneRAtaper(oneRAknots0D2,domain,locs,z,prange, sm, sigeps, knots, r0, r0, ranges, 1,1)
    end
  end

  (loglik_1RBA,T1RBA,TNNZB1)=OneRAblock(oneRAknots0D2,domain, locs,z,prange, sm, sigeps, knots, r0, r0)

###############################################################################
### 1RA slow version
#  for ir0B=1:length(r0Bs),ir0T=1:length(r0Ts),irate= 1:length(rates)

      r0B=0#r0Bs[ir0B]
      r0T=0#r0Ts[ir0T]
      rate=0#rates[irate]
    #(loglik1RBA_S,T1RBA_S,TNNZB1_S)= OneRAblock(domain, locs,z,prange, sm, sigeps, knots, r0B, r0)
    # (loglik_1RTA_S,T1RTA_S,TNNZT1_S) = OneRAtaper(domain, locs,z,prange, sm, sigeps, knots, r0T, r0, ranges, rate, k)
    #try
    #catch e
      #println(e)
    #finally
    #  loglik_1RTA_S=0
    #  T1RTA_S=0
    #  end
    #filePath = "/home/grad/wgong/Spatial/TuningMRTA/parallel/allset0.5_4_"
    fid1 = open(filePath*string(myid())*".csv", "a")
    writecsv(fid1,hcat(r0,J,M,n,sm,prange,k,sigeps,loglikGP0
    ,loglik_MRTA,loglik_1RTA,loglik_MRBA,loglik_1RBA,loglik_1RTA_S,loglik1RBA_S
    ,TGP0,TT,T1RTA,TB,T1RBA,T1RTA_S,T1RBA_S,r0T,rate,r0B,myid()))
    close(fid1)
#  end
end
end

#sim2D(1,1)
#pmap((x1,x2)->sim2D(x1,x2),collect(1:5),repeat([1],inner=[5]))

seed=10#parse(Int64, ARGS[1])+1
#pmap((1,)->sim2D(x1,x2),collect(1:3),repeat([seed],inner=[3]))
sim2D(1,1)#,collect(1:3),repeat([seed],inner=[3])
