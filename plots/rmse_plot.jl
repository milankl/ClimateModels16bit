using FileIO
using Statistics
using PyPlot

var = "u"
R = load("analysis/RMSE_$var.jld2")["RMSE"]

n = size(R)[2]      # number of time steps
nT = size(R)[1]      # number of types
t = (0:(n-1))*0.25  # time array in days, 6h time step

cerr = 0.5      # climatology of the forecast error for long lead time
p = 25

Rp = Array{Float64,3}(undef,2,nT,n)

for i in 1:n
    for j in 1:nT
        Rp[1,j,i] = percentile(R[j,i,:],p)
        Rp[2,j,i] = percentile(R[j,i,:],100-p)
    end
end


ioff()
fig,ax = subplots(1,1,figsize=(6,4))

#ax.plot(t,median(R[1,:,:],dims=2)/cerr,"C0",lw=2,label="persistence")
#ax.plot(t,median(R[2,:,:],dims=2)/cerr,"C0",lw=2,label="discretization error")
#ax.plot(t,median(R[4,:,:],dims=2)/cerr,"C1",lw=2,label="Float32")
ax.plot(t,median(R[4,:,:],dims=2)/cerr,"k",lw=2,label="Float16")
ax.plot(t,median(R[6,:,:],dims=2)/cerr,"#900000",lw=2,label="Posit16_2")
ax.plot(t,median(R[5,:,:],dims=2)/cerr,"C2",lw=2,label="Posit16")
ax.plot(t,median(R[8,:,:],dims=2)/cerr,"C1",lw=2,label="BFloat16/Float32")
ax.plot(t,median(R[7,:,:],dims=2)/cerr,"C0",lw=2,label="Float16/Float32")


ax.fill_between(t,Rp[1,4,:]/cerr,Rp[2,4,:]/cerr,color="k",alpha=0.2)
ax.fill_between(t,Rp[1,5,:]/cerr,Rp[2,5,:]/cerr,color="#50C070",alpha=0.2)
ax.fill_between(t,Rp[1,6,:]/cerr,Rp[2,6,:]/cerr,color="#900000",alpha=0.2)
ax.fill_between(t,Rp[1,7,:]/cerr,Rp[2,7,:]/cerr,color="C0",alpha=0.2)
ax.fill_between(t,Rp[1,8,:]/cerr,Rp[2,8,:]/cerr,color="C1",alpha=0.2)

ax.set_ylabel("normalized RMSE")
ax.set_xlabel("time [days]")
#
ax.set_xlim(0,100)
ax.set_ylim(0,1)
ax.set_title("Forecast error", loc="left")
#
ax.legend()
#
# di = Dict("facecolor"=>"w","alpha"=>0.5,"pad"=>2)
#
# ax[:text](98.5,0.95,"Posit(16,0)",bbox=di,ha="right")
# ax[:text](98.5,0.65,"Float16",bbox=di,ha="right")
# ax[:text](98.5,0.38,"discretisation error",bbox=di,ha="right")
# ax[:text](98.5,0.16,"Posit(16,1)",bbox=di,ha="right")
# ax[:text](98.5,0.03,"Posit(16,2)",bbox=di,ha="right")
# #ax[:text](98.5,0.45,"Posit(16,3)",bbox=di,ha="right")

tight_layout()
savefig("plots/rmse_$var.pdf")
savefig("plots/rmse_$var.png",dpi=300)
close(fig)
