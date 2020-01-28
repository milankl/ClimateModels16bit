using NetCDF
using FileIO
using Printf
using Statistics
using StatsBase

path = "/network/aopp/chaos/pred/kloewer/julsdata/forecast2/"

nn = 8      # number of number types to compare
nt = 404    # number of time steps
ne = 300    # number of ensemble members
nx = 100    # number of grid cells in x
ny = 50     # number of grid cells in y
R = Array{Float64,3}(undef,nn,nt,ne)

vari = "sst"

Δ = 1/40
bins = collect(-Δ/2:Δ:1+Δ/2)

for i in 0:ne-1

    println("Ensemble member $i")

    run_id = "run"*@sprintf("%04d",i)

    # READ TRUTH
    nc = NetCDF.open(joinpath(path,"Float64",run_id,vari*".nc"))
    F64 = nc.vars[vari][:,:,:]
    NetCDF.close(nc)

    # READ F64 with Sadourny/RK3
    nc = NetCDF.open(joinpath(path,"Float64_RK3_Sad",run_id,vari*".nc"))
    F64RK3 = nc.vars[vari][:,:,:]
    # due to differnt time stepping one time step less, copy last time step
    F64RK3 = cat(F64RK3,F64RK3[:,:,end],dims=3)
    NetCDF.close(nc)

    nc = NetCDF.open(joinpath(path,"Float32",run_id,vari*".nc"))
    F32 = nc.vars[vari][:,:,:]
    NetCDF.close(nc)

    nc = NetCDF.open(joinpath(path,"Float16",run_id,vari*".nc"))
    F16 = nc.vars[vari][:,:,:]
    NetCDF.close(nc)

    nc = NetCDF.open(joinpath(path,"Posit16",run_id,vari*".nc"))
    P16 = nc.vars[vari][:,:,:]
    NetCDF.close(nc)

    nc = NetCDF.open(joinpath(path,"Posit16_2",run_id,vari*".nc"))
    P162 = nc.vars[vari][:,:,:]
    NetCDF.close(nc)

    nc = NetCDF.open(joinpath(path,"Float16_32",run_id,vari*".nc"))
    F1632 = nc.vars[vari][:,:,:]
    NetCDF.close(nc)

    nc = NetCDF.open(joinpath(path,"BFloat16_32",run_id,vari*".nc"))
    BF1632 = nc.vars[vari][:,:,:]
    NetCDF.close(nc)

    # Compute - average over x & time (last dimension)
    for (ri,M) in enumerate([F64,F64RK3,F32,F16,P16,P162,F1632,BF1632])
        for it in 1:nt
            #H = fit(Histogram,M[:,:,it][:],bins)
            #R[ri,it,i+1] = entropy(H.weights/sum(H.weights),2)
            R[ri,it,i+1] = entropy(M[:,:,it][:],2)/(nx*ny)*2
        end
    end
end

# OUTPUT
save("entropy.jld2","entropy",R)
