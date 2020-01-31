using NetCDF
using JLD2
using FileIO
using Printf
using Statistics
using StatsBase

pathref = "/network/aopp/chaos/pred/kloewer/julsdata/forecast2/long/"
path = "/network/aopp/chaos/pred/kloewer/julsdata/forecast2/long/"

nn = 5       # number of number types to compare
nt = 1613     # number of time steps
ne = 50     # number of ensemble members
R = Array{Float64,3}(undef,nn,nt,ne)

vari = "eta"

for i in 0:ne-1

    println("Ensemble member $i")

    run_id = "run"*@sprintf("%04d",i)

    # READ TRUTH
    nc = NetCDF.open(joinpath(pathref,"Float64",run_id,vari*".nc"))
    F64 = nc.vars[vari][:,:,:]
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

    # Compute RMSEs - average over space, time is last dim
    R[1,:,i+1] = sqrt.(mean((F64-F16).^2,dims=(1,2)))[1,1,:]
    R[2,:,i+1] = sqrt.(mean((F64-P16).^2,dims=(1,2)))[1,1,:]
    R[3,:,i+1] = sqrt.(mean((F64-P162).^2,dims=(1,2)))[1,1,:]
    R[4,:,i+1] = sqrt.(mean((F64-F1632).^2,dims=(1,2)))[1,1,:]
    R[5,:,i+1] = sqrt.(mean((F64-BF1632).^2,dims=(1,2)))[1,1,:]
end

# OUTPUT
save("RMSE_long_$vari.jld2","RMSE",R)
