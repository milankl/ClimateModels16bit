using NetCDF
using PyPlot
using PyCall
using Statistics
using StatsBase

path = "/network/aopp/chaos/pred/kloewer/julsdata/forecast/tend/"

run_ids = [ "run0000",
            "run0001",
            "run0002",
            "run0003",
            "run0004",
            "run0005",
            "run0006"]

labels = [  "Float64",
            "Float32",
            "BFloat16/Float32",
            "Float16/Float32",
            "Float16",
            "Posit(16,1)",
            "Posit(16,2)"]

colours = [ "C0",
            "#0000A0",
            "grey",
            "C1",
            "k",
            "#50C070",
            "#900000"]

nruns = length(run_ids)

bins = cat(dims=1,[0],10.0.^(-12:0.05:1))
nbins = length(bins)-1

vars = ["u","v","eta","du","dv","deta"]
nvars = 6       # 3 vars + 3 tendencies

H = Array{Float64,3}(undef,nruns,nvars,nbins)
M = Array{Float64,3}(undef,nruns,nvars,3)       # mean and two percentiles
p = 10

for i in 1:nruns
    for j in 1:nvars

        println((i,j))

        # READ DATA
        nc = NetCDF.open(joinpath(path,run_ids[i],vars[j]*".nc"))
        vari = nc.vars[vars[j]][:,:,2:end]
        NetCDF.close(nc)

        absvari = abs.(vari[:])

        hist = fit(Histogram,absvari,bins)
        H[i,j,:] = hist.weights
        M[i,j,1] = mean(absvari)
        M[i,j,2] = percentile(absvari,p)
        M[i,j,3] = percentile(absvari,100-p)
    end
end

## PLOT
ioff()
fig,axs = subplots(2,2,figsize=(7,5),sharex=true,sharey=true)

zp = [1,2,3,4,5,6,7]
z = [4,5,7,6,1,3,2]
zη = [6,7,4,5,1,3,2]

for i in 1:nruns
    axs[1,1].loglog(bins[1:end-1],H[i,1,:],colours[i],label=labels[i],drawstyle="steps-post",zorder=zp[i])
    #axs[1,2].plot(bins[1:end-1],H[i,2,:],colours[i],drawstyle="steps-post",zorder=zp[i])
    axs[1,2].plot(bins[1:end-1],H[i,3,:],colours[i],drawstyle="steps-post",zorder=zp[i])

    axs[2,1].plot(bins[1:end-1],H[i,4,:],colours[i],drawstyle="steps-post",zorder=z[i])
    #axs[2,2].plot(bins[1:end-1],H[i,5,:],colours[i],drawstyle="steps-post",zorder=z[i])
    axs[2,2].plot(bins[1:end-1],H[i,6,:],colours[i],drawstyle="steps-post",zorder=zη[i])
end

# ranges
yscaling(i) = exp((nruns-i+1)/2.3)*6e5

for (a,k) in enumerate([1,3])
    for i in 1:nruns
        axs[1,a].scatter(M[i,k,1],yscaling(i),s=4,c=colours[i])
        axs[1,a].plot(M[i,k,2:3],[1,1]*yscaling(i),colours[i],marker="|",ms=4)

        axs[2,a].scatter(M[i,3+k,1],yscaling(i),s=4,c=colours[i])
        axs[2,a].plot(M[i,3+k,2:3],[1,1]*yscaling(i),colours[i],marker="|",ms=4)
    end

    axs[1,a].text(0.03,-0.02,"//",transform=axs[1,a].transAxes,fontsize=11)
    #axs[1,k].text(0.038,-0.016,"-",transform=axs[1,k].transAxes,color="w",fontweight="bold",fontsize=9)
    axs[2,a].text(0.03,-0.02,"//",transform=axs[2,a].transAxes,fontsize=11)
    #axs[2,k].text(0.038,-0.016,"-",transform=axs[2,k].transAxes,color="w",fontweight="bold",fontsize=9)
end

axs[1,1].set_xlim(bins[2]/2,maximum(bins))
axs[1,1].set_xticks([7.5e-13,1e-10,1e-7,1e-4,1e-1])
axs[1,1].set_xticklabels([0,L"$10^{-10}$",L"$10^{-7}$",L"$10^{-4}$",L"$10^{-1}$"])
axs[1,1].set_ylim(0.5,2e7)

axs[1,1].legend(loc=2,fontsize=8)

axs[1,1].set_title(L"Zonal velocity $u$")
#axs[1,2].set_title(L"Meridional velocity $v$")
axs[1,2].set_title(L"Sea surface height $\eta$")

axs[2,1].set_title(L"Tendency of $u$")
#axs[2,2].set_title(L"Tendency of $v$")
axs[2,2].set_title(L"Tendency of $\eta$")

axs[1,1].set_title("a",loc="right",fontweight="bold")
axs[1,2].set_title("b",loc="right",fontweight="bold")
#axs[1,3].set_title("c",loc="right",fontweight="bold")

axs[2,1].set_title("c",loc="right",fontweight="bold")
axs[2,2].set_title("d",loc="right",fontweight="bold")
#axs[2,3].set_title("f",loc="right",fontweight="bold")

axs[2,1].set_xlabel("value")
axs[2,2].set_xlabel("value")
#axs[2,3].set_xlabel("value")

axs[1,1].set_ylabel("N")
axs[2,1].set_ylabel("N")

tight_layout()
savefig("plots/tendency_hist.png",dpi=300)
savefig("plots/tendency_hist.pdf")
close(fig)
