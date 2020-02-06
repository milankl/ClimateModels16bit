using NetCDF
using PyPlot
using PyCall
using Statistics
using StatsBase

path = "/network/aopp/chaos/pred/kloewer/julsdata/forecast/tendtest/"

run_ids = [ "run0000",
            "run0001",
            "run0002",
            "run0003"]
            # "run0004",
            # "run0005"]

labels = [  "Float64",
            "Float32",
            "BFloat16/Float32",
            "Float16/Float32",
            "Posit16",
            "Posit16_2"]

colours = [ "C0",
            "k",
            "grey",
            "C1",
            "C2",
            "C3"]

nruns = length(run_ids)

bins = 10.0.^(-11:0.05:1)
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
fig,axs = subplots(2,3,figsize=(10,6),sharex=true,sharey=true)

for i in 1:nruns
    axs[1,1].loglog(bins[1:end-1],H[i,1,:],colours[i],label=labels[i],drawstyle="steps-post",zorder=i)
    axs[1,2].plot(bins[1:end-1],H[i,2,:],colours[i],drawstyle="steps-post",zorder=i)
    axs[1,3].plot(bins[1:end-1],H[i,3,:],colours[i],drawstyle="steps-post",zorder=i)

    axs[2,1].plot(bins[1:end-1],H[i,4,:],colours[i],drawstyle="steps-post",zorder=i)
    axs[2,2].plot(bins[1:end-1],H[i,5,:],colours[i],drawstyle="steps-post",zorder=i)
    axs[2,3].plot(bins[1:end-1],H[i,6,:],colours[i],drawstyle="steps-post",zorder=i)
end

# ranges
yscaling(i) = exp(i/1.5)*3e5

for k in 1:3
    for i in 1:nruns
        axs[2,k].scatter(M[i,3+k,1],yscaling(i),s=10,c=colours[i])
        axs[2,k].plot(M[i,3+k,2:3],[1,1]*yscaling(i),colours[i],marker="|")
    end
end

axs[1,1].set_xlim(minimum(bins),maximum(bins))
axs[1,1].set_ylim(0.5,1e7)

axs[1,1].legend(loc=2,fontsize=8)

axs[1,1].set_title(L"Zonal velocity $u$")
axs[1,2].set_title(L"Meridional velocity $v$")
axs[1,3].set_title(L"Sea surface height $\eta$")

axs[2,1].set_title(L"Tendency of $u$")
axs[2,2].set_title(L"Tendency of $v$")
axs[2,3].set_title(L"Tendency of $\eta$")

axs[1,1].set_title("a",loc="right",fontweight="bold")
axs[1,2].set_title("b",loc="right",fontweight="bold")
axs[1,3].set_title("c",loc="right",fontweight="bold")

axs[2,1].set_title("d",loc="right",fontweight="bold")
axs[2,2].set_title("e",loc="right",fontweight="bold")
axs[2,3].set_title("f",loc="right",fontweight="bold")

axs[2,1].set_xlabel("value")
axs[2,2].set_xlabel("value")
axs[2,3].set_xlabel("value")

tight_layout()
savefig("plots/tendency_hist.png",dpi=300)
close(fig)
