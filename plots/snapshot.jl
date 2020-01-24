using NetCDF
using PyPlot
using PyCall

cm = pyimport("cmocean.cm")

path = "/network/aopp/chaos/pred/kloewer/julsdata/forecast2/nx400/"
time_step = 201

run_id1 = "run0001"
run_id2 = "run0002"
run_id3 = "run0003"
run_id4 = "run0004"
run_id5 = "run0001"
run_id6 = "run0001"


# READ DATA
nc = NetCDF.open(joinpath(path,run_id1,"sst.nc"))
sst1 = nc.vars["sst"][:,:,time_step]
x = nc.vars["x"][:]
y = nc.vars["y"][:]
sst1 = reshape(sst1,length(x),length(y))
NetCDF.close(nc)

nc = NetCDF.open(joinpath(path,run_id2,"sst.nc"))
sst2 = reshape(nc.vars["sst"][:,:,time_step],length(x),length(y))
NetCDF.close(nc)

nc = NetCDF.open(joinpath(path,run_id3,"sst.nc"))
sst3 = reshape(nc.vars["sst"][:,:,time_step],length(x),length(y))
NetCDF.close(nc)

nc = NetCDF.open(joinpath(path,run_id4,"sst.nc"))
sst4 = reshape(nc.vars["sst"][:,:,time_step],length(x),length(y))
NetCDF.close(nc)

nc = NetCDF.open(joinpath(path,run_id5,"sst.nc"))
sst5 = reshape(nc.vars["sst"][:,:,time_step],length(x),length(y))
NetCDF.close(nc)

nc = NetCDF.open(joinpath(path,run_id6,"sst.nc"))
sst6 = reshape(nc.vars["sst"][:,:,time_step],length(x),length(y))
NetCDF.close(nc)

## Plot
levs = collect(0:0.05:1)

ioff()
fig,axs = subplots(3,2,figsize=(8,7),sharex=true,sharey=true)

tight_layout(rect=[-.025,0.07,1,0.99],w_pad=0.03,h_pad=0.6)
pos1 = axs[3,1].get_position()
pos2 = axs[3,2].get_position()
cax = fig.add_axes([pos1.x0,0.07,pos2.x1-pos1.x0,0.02])

cmap = cm.tempo_r

q = axs[1,1].contourf(x,y,sst1',levs,cmap=cmap)
axs[2,1].contourf(x,y,sst2',levs,cmap=cmap)
axs[3,1].contourf(x,y,sst3',levs,cmap=cmap)

axs[1,2].contourf(x,y,sst4',levs,cmap=cmap)
axs[2,2].contourf(x,y,sst5',levs,cmap=cmap)
axs[3,2].contourf(x,y,sst6',levs,cmap=cmap)


cb = colorbar(q,cax=cax,orientation="horizontal")
cb.set_label("Tracer concentration")
cb.set_ticks([0,0.2,0.4,0.6,0.8,1])

axs[1,1].set_ylabel(L"$y$")
axs[2,1].set_ylabel(L"$y$")
axs[3,1].set_ylabel(L"$y$")
axs[3,1].set_xlabel(L"$x$")
axs[3,2].set_xlabel(L"$x$")

axs[1,1].set_title("Float64",loc="left",fontsize=8)
axs[2,1].set_title("Posit16",loc="left",fontsize=8)
axs[3,1].set_title("Posit16_2",loc="left",fontsize=8)

axs[1,2].set_title("Float16",loc="left",fontsize=8)
axs[2,2].set_title("Mixed precision: Float16/Float32",loc="left",fontsize=8)
axs[3,2].set_title("Mixed precision: BFloat16/Float32",loc="left",fontsize=8)

axs[1,1].set_title("a",loc="right",fontsize=8,fontweight="bold")
axs[2,1].set_title("b",loc="right",fontsize=8,fontweight="bold")
axs[3,1].set_title("c",loc="right",fontsize=8,fontweight="bold")

axs[1,2].set_title("d",loc="right",fontsize=8,fontweight="bold")
axs[2,2].set_title("e",loc="right",fontsize=8,fontweight="bold")
axs[3,2].set_title("f",loc="right",fontsize=8,fontweight="bold")


axs[1,1].set_xticks([])
axs[1,1].set_yticks([])

savefig("plots/snapshot.png",dpi=300)
close(fig)
