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

vars = ["u","v","eta"]
nvars = 3       # 3 vars + 3 tendencies

nx = 400
ny = 200
nt = 41
tstart = 1
Δ = 5000.0      # km

Uag = Array{Float64,4}(undef,nruns,(nx-2),(ny-2),(nt-tstart+1))
Vag = Array{Float64,4}(undef,nruns,(nx-2),(ny-2),(nt-tstart+1))

g = 10

function coriolis(ϕ::Real,Δy::Real)
    ω = 2π/(24*3600)
    R = 6.371e6
    β = 2*ω/R*cosd(ϕ)
    return 2*ω*sind(ϕ) + β*Δy
end

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

        Ug = -g*f .* Ix(detadx)[:,2:end-1]
        Vg = g*f .* Iy(detady)[2:end-1,:]

        uT = Ix(u[:,:,j])[2:end,2:end-1]
        vT = Iy(v[:,:,j])[2:end-1,:]

        Uag[i,:,:,j] = uT .- Ug
        Vag[i,:,:,j] = vT .- Vg
    end
end

## ageostrophic component projected on the diagonal
UVag = Uag*cosd(-45) .- Vag*sind(-45)

bins = collect(-1.5:0.005:1.5)
nbins = length(bins)-1

H = Array{Float64,2}(undef,nruns,nbins)

for i in 1:nruns
    hist = fit(Histogram,UVag[i,:,:,:][:],bins)
    H[i,:] = hist.weights
end

## PLOT
ioff()
fig,ax = subplots(1,1,figsize=(7,5))

for i in 1:nruns
    ax.plot(bins[1:end-1],H[i,:],colours[i],label=labels[i],drawstyle="steps-post")
end

ax.set_xlim(bins[1],bins[end])

tight_layout()
savefig("plots/ageostrophic.png",dpi=300)
# savefig("plots/ageostrophic.pdf")
close(fig)
