using FileIO
using Statistics
using PyPlot
using PyCall

@pyimport numpy as np

#path = "/network/aopp/cirrus/pred/kloewer/julsdata/forecast2/"
path = "/network/aopp/chaos/pred/kloewer/julsdata/forecast/"
R = load(path*"analysis/RMSE_280.jld2")["RMSE"]

n = size(R)[2]      # number of time steps
t = (0:(n-1))*0.25  # time array in days

cerr = 2      # climatology of the forecast error for long lead time

##
ioff()
fig,ax = subplots(1,1,figsize=(6,4))

ax[:plot](t,median(R[1,:,:],dims=2)/cerr,"C0",lw=2,label="discretization error")
ax[:plot](t,median(R[2,:,:],dims=2)/cerr,"k",lw=2,label="Float16")
#ax[:plot](t,median(R[3,:,:],dims=2)/cerr,"C1",lw=2,label="Float32")
#ax[:plot](t,median(R[4,:,:],dims=2)/cerr,"grey",lw=2,label="persistence")
ax[:plot](t,median(R[5,:,:],dims=2)/cerr,"C1",lw=2,label="Posit(16,0)")
ax[:plot](t,median(R[6,:,:],dims=2)/cerr,"#50C070",lw=2,label="Posit(16,1)")
ax[:plot](t,median(R[7,:,:],dims=2)/cerr,"#700000",lw=2,label="Posit(16,2)")
#ax[:plot](t,median(R[8,:,:],dims=2)/cerr,"C4",lw=2,label="Posit(16,3)")


p = 25

ax[:fill_between](t,np.percentile(R[1,:,:],p,axis=1)/cerr,np.percentile(R[1,:,:],100-p,axis=1)/cerr,facecolor="C0",alpha=0.3)
ax[:fill_between](t,np.percentile(R[2,:,:],p,axis=1)/cerr,np.percentile(R[2,:,:],100-p,axis=1)/cerr,facecolor="k",alpha=0.3)
#ax[:fill_between](t,np.percentile(R[3,:,:],p,axis=1)/cerr,np.percentile(R[3,:,:],100-p,axis=1)/cerr,color="grey",alpha=0.2)
#ax[:fill_between](t,np.percentile(R[4,:,:],p,axis=1)/cerr,np.percentile(R[4,:,:],100-p,axis=1)/cerr,color="grey",alpha=0.2)
ax[:fill_between](t,np.percentile(R[5,:,:],p,axis=1)/cerr,np.percentile(R[5,:,:],100-p,axis=1)/cerr,facecolor="C1",alpha=0.3)
ax[:fill_between](t,np.percentile(R[6,:,:],p,axis=1)/cerr,np.percentile(R[6,:,:],100-p,axis=1)/cerr,facecolor="#50C070",alpha=0.3)
ax[:fill_between](t,np.percentile(R[7,:,:],p,axis=1)/cerr,np.percentile(R[7,:,:],100-p,axis=1)/cerr,facecolor="#700000",alpha=0.2,zorder=-1)
#ax[:fill_between](t,np.percentile(R[7,:,:],p,axis=1)/cerr,np.percentile(R[8,:,:],100-p,axis=1)/cerr,facecolor="C4",alpha=0.3,zorder=-1)

ax[:set_ylabel]("normalized RMSE")
ax[:set_xlabel]("time [days]")

ax[:set_xlim](0,100)
ax[:set_ylim](0,1)
ax[:set_title]("Forecast error", loc="left")

#ax[:legend]()

di = Dict("facecolor"=>"w","alpha"=>0.5,"pad"=>2)

ax[:text](98.5,0.95,"Posit(16,0)",bbox=di,ha="right")
ax[:text](98.5,0.65,"Float16",bbox=di,ha="right")
ax[:text](98.5,0.38,"discretisation error",bbox=di,ha="right")
ax[:text](98.5,0.16,"Posit(16,1)",bbox=di,ha="right")
ax[:text](98.5,0.03,"Posit(16,2)",bbox=di,ha="right")
#ax[:text](98.5,0.45,"Posit(16,3)",bbox=di,ha="right")

tight_layout()
savefig("/home/kloewer/phd/swm_rmse/rmse280_egu.pdf")
close(fig)
