using FileIO
using Statistics
using StatsBase
using PyPlot

Rm = load("analysis/mean_u.jld2")["mean"]
Rv = load("analysis/var_u.jld2")["var"]

nT,ny,ne = size(Rm)      # number of types, grid cells in x, in y, ensemble members

Δx = 20     # grid spacing [km]
x = Δx*collect(0:ny-1).+Δx/2
xf = cat(0,x,1000,dims=1)

p = 25

Rpm = Array{Float64,3}(undef,2,nT,ny)
Rpv = Array{Float64,3}(undef,2,nT,ny)

for i in 1:ny
    for k in 1:nT
        Rpm[1,k,i] = percentile(Rm[k,i,:],p)
        Rpm[2,k,i] = percentile(Rm[k,i,:],100-p)
        Rpv[1,k,i] = percentile(Rv[k,i,:],p)
        Rpv[2,k,i] = percentile(Rv[k,i,:],100-p)
    end
end

u0(x::AbstractArray) = cat(zero(eltype(x)),x,zero(eltype(x)),dims=1)

##
ioff()
fig,(ax1,ax2) = subplots(2,1,figsize=(7,5),sharex=true)

ax1.plot(xf,u0(mean(Rm[1,:,:],dims=2)),"C0",lw=3,label="Float64")
ax1.plot(xf,u0(mean(Rm[3,:,:],dims=2)),"k",lw=1.5,label="Float16")
ax1.plot(xf,u0(mean(Rm[4,:,:],dims=2)),"#50C070",ls="--",lw=1.5,label="Posit(16,1)")
ax1.plot(xf,u0(mean(Rm[5,:,:],dims=2)),"#900000",ls="--",lw=1.5,label="Posit(16,2)")
ax1.plot(xf,u0(mean(Rm[6,:,:],dims=2)),"grey",ls="-.",lw=1.5,label="BFloat16/Float32")
ax1.plot(xf,u0(mean(Rm[7,:,:],dims=2)),"C1",ls="-.",lw=1.5,label="Float16/Float32")
ax1.plot(xf,u0(zero(x)),"k",lw=0.1)

ax1.fill_between(xf,u0(Rpm[1,1,:]),u0(Rpm[2,1,:]),color="C0",alpha=0.1,label="Float64 ensemble")

f64, = ax2.plot(xf,u0(mean(Rv[1,:,:],dims=2)),"C0",lw=2,label="Float64")
f64e = ax2.fill_between(xf,u0(Rpv[1,1,:]),u0(Rpv[2,1,:]),color="C0",alpha=0.1,label="Float64 ensemble")

f16, = ax2.plot(xf,u0(mean(Rv[3,:,:],dims=2)),"k",lw=1.5,label="Float16")
p16, = ax2.plot(xf,u0(mean(Rv[4,:,:],dims=2)),"#50C070",ls="--",lw=1.5,label="Posit(16,1)")
p162, = ax2.plot(xf,u0(mean(Rv[5,:,:],dims=2)),"#900000",ls="--",lw=1.5,label="Posit(16,2)")
bf16, = ax2.plot(xf,u0(mean(Rv[6,:,:],dims=2)),"grey",ls="-.",lw=1.5,label="BFloat16/Float32")
f16m, = ax2.plot(xf,u0(mean(Rv[7,:,:],dims=2)),"C1",ls="-.",lw=1.5,label="Float16/Float32")

labels = ["Float64","Float16","Posit(16,1)","Posit(16,2)","BFloat16/Float32","Float16/Float32"]
ax2.legend([(f64,f64e),f16,p16,p162,bf16,f16m],labels,loc=(0.1,0.02),fontsize=8,ncol=2)

ax1.set_ylabel(L"$u$ [m/s]")
ax2.set_ylabel(L"Variance(u) [$m^2/s^2$]")
ax2.set_xlabel(L"$y$ [km]")

ax1.set_xlim(0,1000)
ax2.set_ylim(0,0.42)
ax2.set_yticks([0,0.1,0.2,0.3,0.4])

ax1.set_title("Mean zonal current", loc="left")
ax2.set_title("Zonal current variability", loc="left")

ax1.set_title("a", fontweight="bold", loc="right")
ax2.set_title("b", fontweight="bold", loc="right")

tight_layout()
savefig("plots/meanvar_u.pdf")
savefig("plots/meanvar_u.png",dpi=300)
close(fig)
