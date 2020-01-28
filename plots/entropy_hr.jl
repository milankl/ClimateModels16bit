using NetCDF
using PyPlot

path = "/network/aopp/chaos/pred/kloewer/julsdata/forecast2/nx400/"

run_id1 = "run0001"
run_id2 = "run0002"
run_id3 = "run0003"
run_id4 = "run0004"
run_id5 = "run0005"
run_id6 = "run0006"

nT = 6
nt = 201
t = (0:(nt-1))*0.25  # time array in days, 6h time step

Δ = 1/100
bins = collect(-Δ/2:Δ:1+Δ/2)

S = Array{Float64,2}(undef,nT,nt)

# READ DATA
nc = NetCDF.open(joinpath(path,run_id1,"sst.nc"))
sst1 = nc.vars["sst"][:,:,:]
NetCDF.close(nc)

nc = NetCDF.open(joinpath(path,run_id2,"sst.nc"))
sst2 = nc.vars["sst"][:,:,:]
NetCDF.close(nc)

nc = NetCDF.open(joinpath(path,run_id3,"sst.nc"))
sst3 = nc.vars["sst"][:,:,:]
NetCDF.close(nc)

nc = NetCDF.open(joinpath(path,run_id4,"sst.nc"))
sst4 = nc.vars["sst"][:,:,:]
NetCDF.close(nc)

nc = NetCDF.open(joinpath(path,run_id5,"sst.nc"))
sst5 = nc.vars["sst"][:,:,:]
NetCDF.close(nc)

nc = NetCDF.open(joinpath(path,run_id6,"sst.nc"))
sst6 = nc.vars["sst"][:,:,:]
NetCDF.close(nc)

for (ri,M) in enumerate([sst1,sst2,sst3,sst4,sst5,sst6])
    for it in 1:nt
        #H = fit(Histogram,M[:,:,it][:],bins)
        #S[ri,it] = entropy(H.weights/sum(H.weights),2)
        S[ri,it] = entropy(M[:,:,it][:],2)/(200*400)*2
        #S[ri,it] = var(M[:,:,it][:])
    end
end


## Plot
ioff()
fig,ax = subplots(1,1,figsize=(6,4))

ax.plot(t,S[1,:,:],"C3",lw=2,label="Float64")
ax.plot(t,S[2,:,:],"k",lw=1,label="Float16")
ax.plot(t,S[3,:,:],"#900000",lw=1,label="Posit16_2")
ax.plot(t,S[4,:,:],"C2",lw=1,label="Posit16")
ax.plot(t,S[5,:,:],"C1",lw=1,label="BFloat16/Float32")
ax.plot(t,S[6,:,:],"C0",lw=1,label="Float16/Float32")

# ax.set_ylabel("normalized RMSE")
ax.set_xlabel("time [days]")
#
ax.set_xlim(0,100)
ax.set_ylim(0,1)
ax.set_title("Mixing entropy", loc="left")
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
savefig("plots/entropy_hr.pdf")
savefig("plots/entropy_hr.png",dpi=300)
close(fig)
