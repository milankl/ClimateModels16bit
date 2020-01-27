using NetCDF
using JLD2
using FileIO
using Printf
using Statistics
using StatsBase

path = "/network/aopp/chaos/pred/kloewer/julsdata/forecast2/"

nn = 8          # number of number types to compare
ny = 50         # number of grid cells in y
ne = 300        # number of ensemble members

Rm = Array{Float64,3}(undef,nn,ny,ne)
Rv = Array{Float64,3}(undef,nn,ny,ne)

nt0 = 200       # take average from timestep nt0

vari = "u"

for i in 0:ne-1

    println("Ensemble member $i")

    run_id = "run"*@sprintf("%04d",i)

    # READ TRUTH
    nc = NetCDF.open(joinpath(path,"Float64",run_id,vari*".nc"))
    F64 = nc.vars[vari][:,:,nt0:end]
    NetCDF.close(nc)

    # Persistence forecast
    F64p = repeat(F64[:,:,1],1,1,size(F64)[3])

    # READ F64 with Sadourny/RK3
    nc = NetCDF.open(joinpath(path,"Float64_RK3_Sad",run_id,vari*".nc"))
    F64RK3 = nc.vars[vari][:,:,nt0:end]
    NetCDF.close(nc)

    nc = NetCDF.open(joinpath(path,"Float32",run_id,vari*".nc"))
    F32 = nc.vars[vari][:,:,nt0:end]
    NetCDF.close(nc)

    nc = NetCDF.open(joinpath(path,"Float16",run_id,vari*".nc"))
    F16 = nc.vars[vari][:,:,nt0:end]
    NetCDF.close(nc)

    nc = NetCDF.open(joinpath(path,"Posit16",run_id,vari*".nc"))
    P16 = nc.vars[vari][:,:,nt0:end]
    NetCDF.close(nc)

    nc = NetCDF.open(joinpath(path,"Posit16_2",run_id,vari*".nc"))
    P162 = nc.vars[vari][:,:,nt0:end]
    NetCDF.close(nc)

    nc = NetCDF.open(joinpath(path,"Float16_32",run_id,vari*".nc"))
    F1632 = nc.vars[vari][:,:,nt0:end]
    NetCDF.close(nc)

    nc = NetCDF.open(joinpath(path,"BFloat16_32",run_id,vari*".nc"))
    BF1632 = nc.vars[vari][:,:,nt0:end]
    NetCDF.close(nc)

    # Compute - average over x & time (last dimension)
    for (ri,M) in enumerate([F64,F64RK3,F32,F16,P16,P162,F1632,BF1632])
        for k in 1:ny
            Rm[ri,k,i+1] = mean(M[:,k,:])
            Rv[ri,k,i+1] = var(M[:,k,:])
        end
    end
end

# OUTPUT
save("mean_$vari.jld2","mean",Rm)
save("var_$vari.jld2","var",Rv)
