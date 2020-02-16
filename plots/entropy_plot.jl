using FileIO
using Statistics
using StatsBase
using PyPlot

R = load("analysis/entropy_LR.jld2")["entropy"]

nT,n,ne = size(R)      # number of types, time steps, ensemble members
p = 25
t = (0:(n-1))*0.5      # time vector

Ridx = BitArray(undef,n)     # index to remove double entries in R, due to SL time stepping
Ridx[1] = true
for i in 2:n
    if R[1,i,1] == R[1,i-1,1]
        Ridx[i] = false
    else
        Ridx[i] = true
    end
end

# low resolution has a different time step, do for that separately
R2 = R[1:end-1,Ridx,:]
t2 = t[Ridx]

Ridx = BitArray(undef,n)     # index to remove double entries in R, due to SL time stepping
Ridx[1] = true
for i in 2:n
    if R[end,i,1] == R[end,i-1,1]
        Ridx[i] = false
    else
        Ridx[i] = true
    end
end

R_LR = R[end,Ridx,:]
t_LR = t[Ridx]
##

nT,n,ne = size(R2)      # number of types, time steps, ensemble members

Rp = Array{Float64,3}(undef,2,nT,n)

for i in 1:n
    for j in 1:nT
        Rp[1,j,i] = percentile(R2[j,i,:],p)
        Rp[2,j,i] = percentile(R2[j,i,:],100-p)
    end
end

n,ne = size(R_LR)      # number of types, time steps, ensemble members
Rp_LR = Array{Float64,2}(undef,2,n)
for i in 1:n
    Rp_LR[1,i] = percentile(R_LR[i,:],p)
    Rp_LR[2,i] = percentile(R_LR[i,:],100-p)
end

##
ioff()
fig,ax = subplots(1,1,figsize=(6,4))

q = 3

ax.plot(t2,mean(R2[1,:,:],dims=2).^q,"C0",lw=3,label="Float64")
ax.plot(t2,mean(R2[2,:,:],dims=2).^q,"k",label="Float16")
ax.plot(t2,mean(R2[3,:,:],dims=2).^q,"#50C070",ls="--",label="Posit(16,1)")
ax.plot(t2,mean(R2[4,:,:],dims=2).^q,"#900000",ls="--",label="Posit(16,2)")
ax.plot(t2,mean(R2[6,:,:],dims=2).^q,"grey",ls="-.",label="BFloat16/Float32")
ax.plot(t2,mean(R2[5,:,:],dims=2).^q,"C1",ls="-.",label="Float16/Float32")
ax.plot(t_LR,mean(R_LR[:,:],dims=2).^q,"C0",ls="--",label="Float64 low res")

ax.fill_between(t2,Rp[1,1,:].^q,Rp[2,1,:].^q,color="C0",alpha=0.2,label="Float64 ensemble",zorder=-1)

ax.set_ylabel("normalized entropy")
ax.set_xlabel("time [days]")

ax.set_xlim(0,300)
ax.set_ylim(0,1)

ytik = [0,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
ax.set_yticks(ytik.^q)
ax.set_yticklabels(string.(ytik))

#ax.set_ylim(0,1)
ax.set_title("Entropy of mixing", loc="left")
#
ax.legend(loc=4,ncol=2)
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
savefig("plots/entropy.pdf")
savefig("plots/entropy.png",dpi=300)
close(fig)
