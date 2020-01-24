using NetCDF
using PyPlot
using PyCall

@pyimport cmocean.cm as cm
#@pyimport matplotlib.pyplot as plt

#plt[:rcParams]["mathtext.fontset"] = "cm"
#plt[:rcParams]["mathtext.rm"] = "serif"

#path = "/network/aopp/chaos/pred/kloewer/forecast2/"
path = "/network/aopp/chaos/pred/kloewer/julsdata/ssthr/"
run_id = "run"*@sprintf("%04d",3)

i = 100

# READ TRUTH
nc = NetCDF.open(path*run_id*"/sst.nc")
sst = reshape(nc.vars["sst"][:,:,i],800,400)
x = nc.vars["x"][:]
y = nc.vars["y"][:]
NetCDF.close(nc)

## SSTA
xx,yy = meshgrid(x,y)
SSTmax = 30.
SSTmin = 0.
SSTw = 1000e3
SSTϕ = 0.5
Ly = 4000e3
sstm = (SSTmin+SSTmax)/2 .+ tanh.(2π*(Ly/(4*SSTw))*(yy/Ly .- SSTϕ))*(SSTmin-SSTmax)/2
ssta = sst - sstm
## Plot

levs = 0:0.2:30
#clevs = 27:.2:30

ioff()
fig,ax = subplots(1,1,figsize=(6,3))
tight_layout(rect=[-0.06,-0.09,0.91,0.99])
pos1 = ax[:get_position]()
cax = fig[:add_axes]([pos1[:x1]+0.01,pos1[:y0],0.02,pos1[:y1]-pos1[:y0]])

q = ax[:contourf](x,y,sst',levs,cmap="magma",extend="both")
#q = ax[:contourf](x,y,sst',levs,colors="k",linewidths=0.2)
cb = colorbar(q,cax=cax)
cb[:set_label]("[°C]")
cb[:set_ticks]([0,5,10,15,20,25,30])
#cb[:set_ticks]([-15,-10,-5,0,5,10,15])

#ax[:set_ylabel](L"$y$")
#ax[:set_xlabel](L"$x$")
ax[:set_xticks]([])
ax[:set_yticks]([])
ax[:set_title]("Temperature anomaly simulated with 64bit floats",loc="left",fontsize=8)

#savefig("/home/kloewer/julia/swm_rmse/snapshot_posit.png",dpi=300)
savefig("/home/kloewer/juls/figs/sst_hr_float64.png",dpi=300)
close(fig)
