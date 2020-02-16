using NetCDF
using FileIO
using Printf
using Statistics
using StatsBase

path = "/network/aopp/chaos/pred/kloewer/julsdata/forecast/"

nn = 1          # number of number types to compare
ny = 25         # number of grid cells in y
ne = 200        # number of ensemble members

Rm = Array{Float64,2}(undef,ny,ne)
Rv = Array{Float64,2}(undef,ny,ne)

nt0 = 100       # take average from timestep nt0

vari = "u"

for i in 0:ne-1

    println("Ensemble member $i")

    run_id = "run"*@sprintf("%04d",i)

    # READ
    nc = NetCDF.open(joinpath(path,"Float64_LR",run_id,vari*".nc"))
    F64LR = nc.vars[vari][:,:,nt0:end]
    NetCDF.close(nc)

    # Compute - average over x & time (last dimension)
    for k in 1:ny
        Rm[k,i+1] = mean(F64LR[:,k,:])
        Rv[k,i+1] = var(F64LR[:,k,:])
    end
end

# OUTPUT
save("mean_LR_$vari.jld2","mean",Rm)
save("var_LR_$vari.jld2","var",Rv)
