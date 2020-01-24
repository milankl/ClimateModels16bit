using NetCDF
using JLD2
using FileIO
using Printf
using Statistics

#path = "/network/aopp/cirrus/pred/kloewer/julsdata/forecast2/"
path = "/network/aopp/chaos/pred/kloewer/julsdata/forecast/"

nn = 8       # number of number types to compare
nt = 404     # number of time steps
ne = 270     # number of ensemble members

R = Array{Float64,3}(undef,nn,nt,ne)

#failed_runs = [9,21,23,24,29,42,43,47,58,68,93,99,109,116,117,122,
                # 139,144,147,150,152,168,190,195,204,211,216,220,223,234,
                # 242,247,251,267,269,271,273,282,283,289,297,301,310,
                # 311,317,319,322,324,328,329,330,345,346,352,356,361,
                # 366,372,376,380,401,406]

#failed_runs = [84,251,252,255,256,257,258,259,260,261,262,263,265,266,267,268,269]
failed_runs = [84,182,183,184,185,186,187,188,192,251,252,255,256,257,258,259,260,261,262,263,265,266,267,268,269]

io = 0  # index for output
ie = 0  # index for runs

while io < ne    # loop over ensemble member
    if ie in failed_runs
        println("Skipping failed ensemble member $ie")
    else
        global io += 1

        println("Ensemble member $ie")

        run_id = "run"*@sprintf("%04d",ie)

        # READ TRUTH
        nc = NetCDF.open(path*"Float64/"*run_id*"/eta.nc")
        F64 = nc.vars["eta"][:,:,:]
        NetCDF.close(nc)

        # Persistence forecast
        F64p = repeat(F64[:,:,1],1,1,size(F64)[3])

        # READ F64 with Sadourny/RK3
        nc = NetCDF.open(path*"Float64Sad/"*run_id*"/eta.nc")
        F64S = nc.vars["eta"][:,:,:]
        #F64S = cat(F64S,F64S[:,:,end],dims=3)   # due to differnt time stepping one time step less, copy...
        NetCDF.close(nc)

        # READ F16
        nc = NetCDF.open(path*"Float16/"*run_id*"/eta.nc")
        F16 = nc.vars["eta"][:,:,:]
        NetCDF.close(nc)

        # READ F32
        nc = NetCDF.open(path*"Float32/"*run_id*"/eta.nc")
        F32 = nc.vars["eta"][:,:,:]
        NetCDF.close(nc)

        # READ POSIT{16,0}
        nc = NetCDF.open(path*"Posit160/"*run_id*"/eta.nc")
        P160 = nc.vars["eta"][:,:,:]
        NetCDF.close(nc)

        # READ POSIT{16,1}
        nc = NetCDF.open(path*"Posit161/"*run_id*"/eta.nc")
        P161 = nc.vars["eta"][:,:,:]
        NetCDF.close(nc)

        # READ POSIT{16,2}
        nc = NetCDF.open(path*"Posit162/"*run_id*"/eta.nc")
        P162 = nc.vars["eta"][:,:,:]
        NetCDF.close(nc)

        # READ POSIT{16,3}
        nc = NetCDF.open(path*"Posit163/"*run_id*"/eta.nc")
        P163 = nc.vars["eta"][:,:,:]
        NetCDF.close(nc)

        # Compute RMSEs - average over space, time is last dim
        R[1,:,io] = sqrt.(mean((F64-F64S).^2,dims=(1,2)))[1,1,:]
        R[2,:,io] = sqrt.(mean((F64-F16).^2,dims=(1,2)))[1,1,:]
        R[3,:,io] = sqrt.(mean((F64-F32).^2,dims=(1,2)))[1,1,:]
        R[4,:,io] = sqrt.(mean((F64-F64p).^2,dims=(1,2)))[1,1,:]
        R[5,:,io] = sqrt.(mean((F64-P160).^2,dims=(1,2)))[1,1,:]
        R[6,:,io] = sqrt.(mean((F64-P161).^2,dims=(1,2)))[1,1,:]
        R[7,:,io] = sqrt.(mean((F64-P162).^2,dims=(1,2)))[1,1,:]
        R[8,:,io] = sqrt.(mean((F64-P163).^2,dims=(1,2)))[1,1,:]
    end

    global ie += 1

end

# OUTPUT
save(path*"analysis/RMSE_270P3.jld2","RMSE",R)
