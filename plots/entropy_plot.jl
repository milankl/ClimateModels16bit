using FileIO
using Statistics
using PyPlot

R = load("analysis/entropy.jld2")["entropy"]*2

nT,n,ne = size(R)      # number of types, time steps, ensemble members
p = 25
t = (0:(n-1))*0.25      # time vector

Ridx = BitArray(undef,n)     # index to remove double entries in R, due to SL time stepping
Ridx[1] = true
for i in 2:n
    if R[1,i,1] == R[1,i-1,1]
        Ridx[i] = false
    else
        Ridx[i] = true
    end
end

R = R[:,Ridx,:]
t = t[Ridx]
nT,n,ne = size(R)      # number of types, time steps, ensemble members

Rp = Array{Float64,3}(undef,2,nT,n)

for i in 1:n
    for j in 1:nT
        Rp[1,j,i] = percentile(R[j,i,:],p)
        Rp[2,j,i] = percentile(R[j,i,:],100-p)
    end
end

##
ioff()
fig,ax = subplots(1,1,figsize=(6,4))

q = 2

ax.plot(t,median(R[1,:,:],dims=2).^q,"C3",lw=2,label="Float64",zorder=4)
ax.plot(t,median(R[4,:,:],dims=2).^q,"k",label="Float16")
ax.plot(t,median(R[5,:,:],dims=2).^q,"C2",label="Posit(16,1)")
ax.plot(t,median(R[6,:,:],dims=2).^q,"#900000",label="Posit(16,2)")
ax.plot(t,median(R[8,:,:],dims=2).^q,"C1",label="BFloat16/Float32")
ax.plot(t,median(R[7,:,:],dims=2).^q,"C0",label="Float16/Float32")

ax.fill_between(t,Rp[1,1,:].^q,Rp[2,1,:].^q,color="C3",alpha=0.2,label="Float64 ensemble",zorder=-1)

ax.set_ylabel("normalized entropy")
ax.set_xlabel("time [days]")

ax.set_xlim(0,100)

ytik = [0,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
ax.set_yticks(ytik.^q)
ax.set_yticklabels(string.(ytik))

#ax.set_ylim(0,1)
ax.set_title("Mixing entropy", loc="left")
#
ax.legend(loc=2,ncol=2)
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