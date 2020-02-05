using FileIO
using Statistics
using StatsBase
using PyPlot

var = "eta"
R = load("analysis/RMSE_long_$var.jld2")["RMSE"]
R2 = load("analysis/RMSE_comm_$var.jld2")["RMSE"]

n = size(R)[2]      # number of time steps
nT = size(R)[1]      # number of types
t = (0:(n-1))*0.25  # time array in days, 6h time step

n2 = size(R2)[2]      # number of time steps
nT2 = size(R2)[1]      # number of types
t2 = (0:(n2-1))*0.5  # time array in days, 12h time step

cerr = 1.2      # climatology of the forecast error for long lead time
p = 25

Rp = Array{Float64,3}(undef,2,nT,n)
Rp2 = Array{Float64,3}(undef,2,nT2,n2)

for i in 1:n
    for j in 1:nT
        Rp[1,j,i] = percentile(R[j,i,:],p)
        Rp[2,j,i] = percentile(R[j,i,:],100-p)
    end
end

for i in 1:n2
    for j in 1:nT2
        Rp2[1,j,i] = percentile(R2[j,i,:],p)
        Rp2[2,j,i] = percentile(R2[j,i,:],100-p)
    end
end
##
ioff()
fig,(ax1,ax2) = subplots(1,2,figsize=(8,4),sharey=true,sharex=true)

ax1.semilogy(t,median(R[1,:,:],dims=2)/cerr,"k",label="Float16")
ax1.semilogy(t,median(R[2,:,:],dims=2)/cerr,"#50C070",ls="--",label="Posit(16,1)")
ax1.semilogy(t,median(R[3,:,:],dims=2)/cerr,"#900000",ls="--",label="Posit(16,2)")
ax1.semilogy(t,median(R[5,:,:],dims=2)/cerr,"grey",ls="-.",label="BFloat16,Float32 mixed")
ax1.semilogy(t,median(R[4,:,:],dims=2)/cerr,"C1",ls="-.",label="Float16,Float32 mixed")

ax1.fill_between(t,Rp[1,1,:]/cerr,Rp[2,1,:]/cerr,color="k",alpha=0.2)
ax1.fill_between(t,Rp[1,2,:]/cerr,Rp[2,2,:]/cerr,color="#50C070",alpha=0.2)
ax1.fill_between(t,Rp[1,3,:]/cerr,Rp[2,3,:]/cerr,color="#900000",alpha=0.2)
ax1.fill_between(t,Rp[1,5,:]/cerr,Rp[2,5,:]/cerr,color="grey",alpha=0.2)
ax1.fill_between(t,Rp[1,4,:]/cerr,Rp[2,4,:]/cerr,color="C1",alpha=0.2)

ax1.set_ylabel("normalized RMSE")
ax1.set_xlabel("time [days]")
ax2.set_xlabel("time [days]")

ax1.set_xlim(0,250)
ax1.set_ylim(0,1)
ax1.set_title("Forecast error:\n16bit calculations", loc="left")
ax1.set_title("a", loc="right", fontweight="bold")

ax1.legend(loc=4,fontsize=8)

ax2.semilogy(t2,median(R2[1,:,:],dims=2)/cerr,"k",label="Float16")
ax2.semilogy(t2,median(R2[2,:,:],dims=2)/cerr,"grey",ls="-",label="BFloat16")
ax2.semilogy(t2,median(R2[3,:,:],dims=2)/cerr,"#50C070",ls="--",label="Posit(16,1)")
ax2.semilogy(t2,median(R2[4,:,:],dims=2)/cerr,"#900000",ls="--",label="Posit(16,2)")
ax2.semilogy(t2,median(R2[5,:,:],dims=2)/cerr,"C4",ls="--",label="Posit(8,0)")

ax2.fill_between(t2,Rp2[1,1,:]/cerr,Rp2[2,1,:]/cerr,color="k",alpha=0.2)
ax2.fill_between(t2,Rp2[1,2,:]/cerr,Rp2[2,2,:]/cerr,color="grey",alpha=0.2)
ax2.fill_between(t2,Rp2[1,3,:]/cerr,Rp2[2,3,:]/cerr,color="#50C070",alpha=0.2)
ax2.fill_between(t2,Rp2[1,4,:]/cerr,Rp2[2,4,:]/cerr,color="#900000",alpha=0.2)
ax2.fill_between(t2,Rp2[1,5,:]/cerr,Rp2[2,5,:]/cerr,color="C4",alpha=0.2)

ax2.set_title("Forecast error:\n16 or 8bit communication", loc="left")
ax2.set_title("b", loc="right", fontweight="bold")

ax2.legend(loc=4,fontsize=8)

tight_layout()
savefig("plots/rmse_$var.pdf")
savefig("plots/rmse_$var.png",dpi=300)
close(fig)
