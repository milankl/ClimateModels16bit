using Lorenz63
using PyPlot
using SoftPosit
using BFloat16s

N = 100_000

xyz_ini = [7.88,5.39,20.67]

Δt = 0.008

XYZ_F64 = L63(Float64,N=N,xyz=xyz_ini,Δt=Δt)
XYZ_F16 = L63(Float16,N=N,xyz=xyz_ini,Δt=Δt)
XYZ_P16 = L63(Posit16,N=N,xyz=xyz_ini,Δt=Δt)
XYZ_BF16 = L63(BFloat16,N=N,xyz=xyz_ini,Δt=Δt)

# PLOT
ioff()
fig,axs = subplots(2,2,sharex=true,sharey=true)

axs[1,1].plot(XYZ_F64[1,:],XYZ_F64[3,:],"k",lw=0.2)
axs[1,2].plot(XYZ_F16[1,:],XYZ_F16[3,:],"k",lw=0.2)
axs[2,1].plot(XYZ_P16[1,:],XYZ_P16[3,:],"k",lw=0.2)
axs[2,2].plot(XYZ_BF16[1,:],XYZ_BF16[3,:],"k",lw=0.2)

axs[1,1].set_title("Float64",loc="left")
axs[1,2].set_title("Float16",loc="left")
axs[2,1].set_title("Posit16",loc="left")
axs[2,2].set_title("BFloat16",loc="left")

axs[1,1].set_title("a",loc="right",fontweight="bold")
axs[1,2].set_title("b",loc="right",fontweight="bold")
axs[2,1].set_title("c",loc="right",fontweight="bold")
axs[2,2].set_title("d",loc="right",fontweight="bold")

axs[1,1].set_xticks([])
axs[1,1].set_yticks([])

axs[2,1].set_xlabel(L"$x$")
axs[2,2].set_xlabel(L"$x$")
axs[2,1].set_ylabel(L"$z$")
axs[1,1].set_ylabel(L"$z$")

tight_layout()
savefig("plots/lorenz_attractor.png",dpi=300)
close(fig)
