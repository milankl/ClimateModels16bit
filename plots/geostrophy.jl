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
            "run0006"]

labels = [  "Float64",
            "Float16/Float32",
            "Float16",
            "Posit(16,1)",
            "Posit(16,2)",
            "Float32",
            "BFloat16/Float32"]

colours = [ "C0",
            "C1",
            "k",
            "#50C070",
            "#900000",
            "#0000A0",
            "grey"]

nruns = length(run_ids)

vars = ["u","v","eta"]
nvars = 3       # 3 vars + 3 tendencies

nx = 400
ny = 200
nt = 201
tstart = 100
Δ = 5000.0      # km

Ug = Array{Float64,4}(undef,nruns,(nx-2),(ny-2),(nt-tstart+1))
Vg = Array{Float64,4}(undef,nruns,(nx-2),(ny-2),(nt-tstart+1))

Ugpar = Array{Float64,4}(undef,nruns,(nx-2),(ny-2),(nt-tstart+1))
Ugper = Array{Float64,4}(undef,nruns,(nx-2),(ny-2),(nt-tstart+1))

Uag = Array{Float64,4}(undef,nruns,(nx-2),(ny-2),(nt-tstart+1))
Vag = Array{Float64,4}(undef,nruns,(nx-2),(ny-2),(nt-tstart+1))

Uagpar = Array{Float64,4}(undef,nruns,(nx-2),(ny-2),(nt-tstart+1))
Uagper = Array{Float64,4}(undef,nruns,(nx-2),(ny-2),(nt-tstart+1))

g = 10

function coriolis(ϕ::Real,Δy::Real)
    ω = 2π/(24*3600)
    R = 6.371e6
    β = 2*ω/R*cosd(ϕ)
    return 2*ω*sind(ϕ) + β*Δy
end

R(ϕ) = [cosd(ϕ) -sind(ϕ); sind(ϕ) cosd(ϕ)]
R90 = R(90)

y = collect((-ny/2+1.5)*Δ:Δ:(ny/2-1.5)*Δ)
f = Array{Float64,2}(undef,nx-2,ny-2)
f0 = coriolis(45,0)

for (i,iy) in enumerate(y)
    f[:,i] .= coriolis(45,iy)
end

for i in 1:nruns

    nc = NetCDF.open(joinpath(path,run_ids[i],"u.nc"))
    u = nc.vars["u"][:,:,tstart:end]
    NetCDF.close(nc)

    nc = NetCDF.open(joinpath(path,run_ids[i],"v.nc"))
    v = nc.vars["v"][:,:,tstart:end]
    NetCDF.close(nc)

    nc = NetCDF.open(joinpath(path,run_ids[i],"eta.nc"))
    eta = nc.vars["eta"][:,:,tstart:end]
    NetCDF.close(nc)

    for j in 1:(nt-tstart+1)

        detadx = ∂x(eta[:,:,j],Δ)
        detady = ∂y(eta[:,:,j],Δ)

        ug = -g./f .* Iy(detady)[2:end-1,:]
        vg = g./f .* Ix(detadx)[:,2:end-1]

        uT = Ix(u[:,:,j])[2:end,2:end-1]
        vT = Iy(v[:,:,j])[2:end-1,:]

        Uag[i,:,:,j] = uT .- ug
        Vag[i,:,:,j] = vT .- vg

        # Rotate the vectors in a flow parallel and a flow perpendicular component
        # uag_par = (uag'*u)/(u'*u) * u
        # uag_per = (uag'*u90)/(u90'*u90) * u90

        for l in 1:nx-2
            for n in 1:ny-2
                uv = [uT[l,n],vT[l,n]]
                uv90 = R90*uv
                uvag = [Uag[i,l,n,j],Vag[i,l,n,j]]

                Uagpar[i,l,n,j] = (uvag'*uv)/sqrt(uv'*uv)
                Uagper[i,l,n,j] = (uvag'*uv90)/sqrt(uv90'*uv90)
            end
        end

        for l in 1:nx-2
            for n in 1:ny-2
                uv = [uT[l,n],vT[l,n]]
                uv90 = R90*uv
                uvg = [ug[l,n],vg[l,n]]

                Ugpar[i,l,n,j] = (uvg'*uv)/sqrt(uv'*uv)
                Ugper[i,l,n,j] = (uvg'*uv90)/sqrt(uv90'*uv90)
            end
        end

    end
end

## ageostrophic component projected on the diagonal
bins1 = collect(0:0.004:2.4)
nbins1 = length(bins1)-1

Hgpar = Array{Float64,2}(undef,nruns,nbins1)
Hgper = Array{Float64,2}(undef,nruns,nbins1)

bins2 = collect(-.3:0.001:.3)
nbins2 = length(bins2)-1

Hagpar = Array{Float64,2}(undef,nruns,nbins2)
Hagper = Array{Float64,2}(undef,nruns,nbins2)


for i in 1:nruns
    hist = fit(Histogram,Ugpar[i,:,:,:][:],bins1)
    Hgpar[i,:] = hist.weights

    hist = fit(Histogram,Ugper[i,:,:,:][:],bins1)
    Hgper[i,:] = hist.weights

    hist = fit(Histogram,Uagpar[i,:,:,:][:],bins2)
    Hagpar[i,:] = hist.weights

    hist = fit(Histogram,Uagper[i,:,:,:][:],bins2)
    Hagper[i,:] = hist.weights
end

## PLOT
ioff()
fig,(ax1,ax2) = subplots(2,1,figsize=(6,5))

for i in [1,6,3,4,5,2,7]
    ax1.plot(bins1[1:end-1],Hgpar[i,:]/1e4,colours[i],label=labels[i],drawstyle="steps-post")
    ax2.plot(bins2[1:end-1],Hpar[i,:]/1e5,colours[i],label=labels[i],drawstyle="steps-post")
end


ax1.set_xlim(bins1[1],bins1[end])
ax2.set_xlim(bins2[1],bins2[end])
ax1.set_ylabel(L"N [10$^4$]")
ax2.set_ylabel(L"N [10$^5$]")
# ax1.set_xlabel("[m/s]")
ax2.set_xlabel("[m/s]")

ax1.set_yticks([0,1,2,3,4])
ax2.set_yticks([0,1,2])

ax1.set_title("Geostrophic velocity, flow parallel component", loc="left")
ax2.set_title("Ageostrophic velocity, flow parallel component", loc="left")
ax1.set_title("a", loc="right", fontweight="bold")
ax2.set_title("b", loc="right", fontweight="bold")

ax2.legend(loc=2,fontsize=8)

tight_layout()
savefig("plots/ageostrophic.png",dpi=300)
# savefig("plots/ageostrophic.pdf")
close(fig)
