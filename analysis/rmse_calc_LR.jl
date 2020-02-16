using NetCDF
using JLD2
using FileIO
using Printf
using Statistics
using StatsBase
using Interpolations

path = "/network/aopp/chaos/pred/kloewer/julsdata/forecast/"

nn = 8       # number of number types to compare
nt = 601     # number of time steps
ne = 200     # number of ensemble members
R = Array{Float64,3}(undef,nn,nt,ne)

vari = "eta"

for i in 0:ne-1

    println("Ensemble member $i")

    run_id = "run"*@sprintf("%04d",i)

    # READ TRUTH
    nc = NetCDF.open(joinpath(path,"Float64",run_id,vari*".nc"))
    F64 = nc.vars[vari][:,:,:]
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

    nc = NetCDF.open(joinpath(path,"Posit32",run_id,vari*".nc"))
    P32 = nc.vars[vari][:,:,:]
    NetCDF.close(nc)

    nc = NetCDF.open(joinpath(path,"Float64_LR",run_id,vari*".nc"))
    F64LR = nc.vars[vari][:,:,1:nt]
    NetCDF.close(nc)

    # Compute RMSEs - average over space, time is last dim
    R[1,:,i+1] = sqrt.(mean((F64-F32).^2,dims=(1,2)))[1,1,:]
    R[2,:,i+1] = sqrt.(mean((F64-F16).^2,dims=(1,2)))[1,1,:]
    R[3,:,i+1] = sqrt.(mean((F64-P16).^2,dims=(1,2)))[1,1,:]
    R[4,:,i+1] = sqrt.(mean((F64-P162).^2,dims=(1,2)))[1,1,:]
    R[5,:,i+1] = sqrt.(mean((F64-F1632).^2,dims=(1,2)))[1,1,:]
    R[6,:,i+1] = sqrt.(mean((F64-BF1632).^2,dims=(1,2)))[1,1,:]
    R[7,:,i+1] = sqrt.(mean((F64-P32).^2,dims=(1,2)))[1,1,:]

    # interpolate F64 on F64LR
    nx,ny = size(F64)
    x_T = collect(0.5:nx-0.5)
    y_T = collect(0.5:ny-0.5)

    # new grids half the resolution
    Δx,Δy = 2,2

    x_T_new = collect(Δx/2:Δx:nx-Δx/2)
    y_T_new = collect(Δy/2:Δy:ny-Δy/2)

    for it in 1:nt
        η = F64[:,:,it]
        η_itp = interpolate((x_T,y_T),η,Gridded(Linear()))
        F64i = η_itp(x_T_new,y_T_new)
        R[8,it,i+1] = sqrt(mean((F64i-F64LR[:,:,it]).^2))
    end

end

# OUTPUT
save("RMSE_LR_$vari.jld2","RMSE",R)
