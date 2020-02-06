using NetCDF
using FileIO
using Printf
using Statistics
using StatsBase

path = "/network/aopp/chaos/pred/kloewer/julsdata/forecast/comm/"

nn = 5          # number of number types to compare
ny = 50         # number of grid cells in y
ne = 200        # number of ensemble members

Rm = Array{Float64,3}(undef,nn,ny,ne)
Rv = Array{Float64,3}(undef,nn,ny,ne)

nt0 = 100       # take average from timestep nt0

vari = "u"

for i in 0:ne-1

    println("Ensemble member $i")

    run_id = "run"*@sprintf("%04d",i)

    # READ TRUTH
    nc = NetCDF.open(joinpath(path,"Float16",run_id,vari*".nc"))
    F16 = nc.vars[vari][:,:,nt0:end]
    NetCDF.close(nc)

    nc = NetCDF.open(joinpath(path,"Posit16",run_id,vari*".nc"))
    P16 = nc.vars[vari][:,:,nt0:end]
    NetCDF.close(nc)

    nc = NetCDF.open(joinpath(path,"Posit16_2",run_id,vari*".nc"))
    P162 = nc.vars[vari][:,:,nt0:end]
    NetCDF.close(nc)

    nc = NetCDF.open(joinpath(path,"BFloat16",run_id,vari*".nc"))
    BF16 = nc.vars[vari][:,:,nt0:end]
    NetCDF.close(nc)

    nc = NetCDF.open(joinpath(path,"Posit8",run_id,vari*".nc"))
    P8 = nc.vars[vari][:,:,nt0:end]
    NetCDF.close(nc)

    # Compute - average over x & time (last dimension)
    for (ri,M) in enumerate([F16,P16,P162,BF16,P8])
        for k in 1:ny
            Rm[ri,k,i+1] = mean(M[:,k,:])
            Rv[ri,k,i+1] = var(M[:,k,:])
        end
    end
end

# OUTPUT
save("mean_comm_$vari.jld2","mean",Rm)
save("var_comm_$vari.jld2","var",Rv)
