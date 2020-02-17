using NetCDF
using PyPlot
using PyCall
using Statistics
using StatsBase
using ShallowWaters

path = "/network/aopp/chaos/pred/kloewer/julsdata/forecast/geostrophic/"

run_ids = [ "run0000",
            "run0001",
            "run0002",
            "run0003",
            "run0004",
            "run0005",
            "run0006",
            "run0007"]

labels = [  "Float64",
            "Float16/Float32",
            "Float16",
            "Posit(16,1)",
            "Posit(16,2)",
            "Float32",
            "BFloat16/Float32",
            "Float64 low resolution"]

colours = [ "C0",
            "C1",
            "k",
            "#50C070",
            "#900000",
            "#0000A0",
            "grey",
            "#E0E020"]

nruns = length(run_ids)

nx = 400
ny = 200
nt = 201
tstart = 20

Ugpar = Array{Float64,4}(undef,nruns-1,(nx-2),(ny-2),(nt-tstart+1))
Ugper = Array{Float64,4}(undef,nruns-1,(nx-2),(ny-2),(nt-tstart+1))

Uagpar = Array{Float64,4}(undef,nruns-1,(nx-2),(ny-2),(nt-tstart+1))
Uagper = Array{Float64,4}(undef,nruns-1,(nx-2),(ny-2),(nt-tstart+1))

UgparLR = Array{Float64,3}(undef,(200-2),(100-2),(nt-tstart+1))
UgperLR = Array{Float64,3}(undef,(200-2),(100-2),(nt-tstart+1))

UagparLR = Array{Float64,3}(undef,(200-2),(100-2),(nt-tstart+1))
UagperLR = Array{Float64,3}(undef,(200-2),(100-2),(nt-tstart+1))

g = 10

function coriolis(ϕ::Real,Δy::Real)
    ω = 2π/(24*3600)
    R = 6.371e6
    β = 2*ω/R*cosd(ϕ)
    return 2*ω*sind(ϕ) + β*Δy
end

R(ϕ) = [cosd(ϕ) -sind(ϕ); sind(ϕ) cosd(ϕ)]
R90 = R(90)

for i in 1:nruns

    nc = NetCDF.open(joinpath(path,run_ids[i],"u.nc"))
    u = nc.vars["u"][:,:,tstart:nt]
    NetCDF.close(nc)

    nc = NetCDF.open(joinpath(path,run_ids[i],"v.nc"))
    v = nc.vars["v"][:,:,tstart:nt]
    NetCDF.close(nc)

    nc = NetCDF.open(joinpath(path,run_ids[i],"eta.nc"))
    eta = nc.vars["eta"][:,:,tstart:nt]
    NetCDF.close(nc)

    nx,ny,ntime = size(eta)

    Δ = 2000e3/nx

    y = collect((-ny/2+1.5)*Δ:Δ:(ny/2-1.5)*Δ)
    f = Array{Float64,2}(undef,nx-2,ny-2)
    f0 = coriolis(45,0)

    for (ii,iiy) in enumerate(y)
        f[:,ii] .= coriolis(45,iiy)
    end

    for j in 1:ntime

        detadx = ∂x(eta[:,:,j],Δ)
        detady = ∂y(eta[:,:,j],Δ)

        ug = -g./f .* Iy(detady)[2:end-1,:]
        vg = g./f .* Ix(detadx)[:,2:end-1]

        uT = Ix(u[:,:,j])[2:end,2:end-1]
        vT = Iy(v[:,:,j])[2:end-1,:]

        uag = uT .- ug
        vag = vT .- vg

        # Rotate the vectors in a flow parallel and a flow perpendicular component
        # uag_par = (uag'*u)/(u'*u) * u
        # uag_per = (uag'*u90)/(u90'*u90) * u90

        # AGEOSTROPHIC COMPONENT
        for l in 1:nx-2
            for n in 1:ny-2
                uv = [uT[l,n],vT[l,n]]
                uv90 = R90*uv
                uvag = [uag[l,n],vag[l,n]]

                if run_ids[i] == "run0007"  # smaller resolution requires differnt storing matrix
                    UagparLR[l,n,j] = (uvag'*uv)/sqrt(uv'*uv)
                    UagperLR[l,n,j] = (uvag'*uv90)/sqrt(uv90'*uv90)
                else
                    Uagpar[i,l,n,j] = (uvag'*uv)/sqrt(uv'*uv)
                    Uagper[i,l,n,j] = (uvag'*uv90)/sqrt(uv90'*uv90)
                end
            end
        end

        # GEOSTROPHIC COMPONENT
        for l in 1:nx-2
            for n in 1:ny-2
                uv = [uT[l,n],vT[l,n]]
                uv90 = R90*uv
                uvg = [ug[l,n],vg[l,n]]

                if run_ids[i] == "run0007"  # smaller resolution requires differnt storing matrix
                    UgparLR[l,n,j] = (uvg'*uv)/sqrt(uv'*uv)
                    UgperLR[l,n,j] = (uvg'*uv90)/sqrt(uv90'*uv90)
                else
                    Ugpar[i,l,n,j] = (uvg'*uv)/sqrt(uv'*uv)
                    Ugper[i,l,n,j] = (uvg'*uv90)/sqrt(uv90'*uv90)
                end
            end
        end

    end
end

## ageostrophic component projected on the diagonal
bins1 = collect(-0.2:0.005:2.8)
nbins1 = length(bins1)-1

Hgpar = Array{Float64,2}(undef,nruns,nbins1)
Hgper = Array{Float64,2}(undef,nruns,nbins1)
Mgpar = Array{Float64,2}(undef,nruns,3)
Mgper = Array{Float64,2}(undef,nruns,3)

bins2 = collect(-.3:0.001:.3)
nbins2 = length(bins2)-1

Hagpar = Array{Float64,2}(undef,nruns,nbins2)
Hagper = Array{Float64,2}(undef,nruns,nbins2)
Magpar = Array{Float64,2}(undef,nruns,3)
Magper = Array{Float64,2}(undef,nruns,3)

p = 10

for i in 1:nruns

    if i != 8
        geo_parallel = Ugpar[i,:,:,:][:]
        geo_perpendi = Ugper[i,:,:,:][:]
        ageo_parallel = Uagpar[i,:,:,:][:]
        ageo_perpendi = Uagper[i,:,:,:][:]
        f = 1   # histogram factor to account for different grid sizes
    else
        geo_parallel = UgparLR[:,:,:][:]
        geo_perpendi = UgperLR[:,:,:][:]
        ageo_parallel = UagparLR[:,:,:][:]
        ageo_perpendi = UagperLR[:,:,:][:]
        f = 4
    end

    hist = fit(Histogram,geo_parallel,bins1)
    Hgpar[i,:] = hist.weights*f

    hist = fit(Histogram,geo_perpendi,bins1)
    Hgper[i,:] = hist.weights*f

    hist = fit(Histogram,ageo_parallel,bins2)
    Hagpar[i,:] = hist.weights*f

    hist = fit(Histogram,ageo_perpendi,bins2)
    Hagper[i,:] = hist.weights*f

    Mgpar[i,1] = mean(geo_parallel)
    Mgper[i,1] = mean(geo_perpendi)
    Magpar[i,1] = mean(ageo_parallel)
    Magper[i,1] = mean(ageo_perpendi)

    Mgpar[i,2] = percentile(geo_parallel,p)
    Mgper[i,2] = percentile(geo_perpendi,p)
    Magpar[i,2] = percentile(ageo_parallel,p)
    Magper[i,2] = percentile(ageo_perpendi,p)

    Mgpar[i,3] = percentile(geo_parallel,100-p)
    Mgper[i,3] = percentile(geo_perpendi,100-p)
    Magpar[i,3] = percentile(ageo_parallel,100-p)
    Magper[i,3] = percentile(ageo_perpendi,100-p)
end

## PLOT

lws = [3,1.5,1.5,1.5,1.5,1.5,1.5,1]

ioff()
fig,(ax1,ax2) = subplots(1,2,figsize=(8,4))

for (Δ,i) in enumerate([1,6,3,4,5,2,7,8])
    ax1.plot(bins1[1:end-1],Hgpar[i,:]/1e4,colours[i],label=labels[i],lw=lws[i],drawstyle="steps-post")
    ax2.plot(bins2[1:end-1],Hagpar[i,:]/1e5,colours[i],label=labels[i],lw=lws[i],drawstyle="steps-post")

    ax1.scatter(Mgpar[i,1],10-Δ/6,s=4,c=colours[i])
    ax1.plot(Mgpar[i,2:3],[1,1]*(10-Δ/6),colours[i],marker="|",ms=4)

    ax2.scatter(Magpar[i,1],4.3-Δ/12,s=4,c=colours[i])
    ax2.plot(Magpar[i,2:3],[1,1]*(4.3-Δ/12),colours[i],marker="|",ms=4)
end

ax1.set_xlim(bins1[1],bins1[end])
ax2.set_xlim(bins2[1],bins2[end])
ax1.set_ylabel(L"N [10$^4$]")
ax2.set_ylabel(L"N [10$^5$]")
ax1.set_xlabel("[m/s]")
ax2.set_xlabel("[m/s]")

ax1.set_ylim(-0.25,10)
ax2.set_ylim(-0.125,4.3)

ax1.set_title("Geostrophic velocity", loc="left")
ax2.set_title("Ageostrophic velocity", loc="left")
ax1.set_title("a", loc="right", fontweight="bold")
ax2.set_title("b", loc="right", fontweight="bold")

ax1.legend(loc=(0.47,0.4),fontsize=8)

tight_layout()
savefig("plots/ageostrophic.png",dpi=300)
# savefig("plots/ageostrophic.pdf")
close(fig)
