using Lorenz63

function lorenz_data(   N::Int,
                        s::Real=1.0,
                        σ::Real=10.,
                        ρ::Real=28.0,
                        β::Real=8/3,
                        Δt::Real=0.002)

    xyz = L63(N=N-1,σ=σ,ρ=ρ,β=β,Δt=Δt)
    x,y,z = s*xyz[1,:], s*xyz[2,:], s*xyz[3,:]
    RKs = 1/6.

    # create different terms, with scaling
    s_inv = 1/s

    dx = y .- x
    zs = z*s_inv
    ρz = ρ .- zs
    xρz = x .* ρz
    dy = xρz .- y
    ys = y*s_inv
    xy = x .* ys
    βz = β * z
    dz = xy .- βz

    # σ*Δt*RKs,Δt*RKs are precomputed - do here in one step
    tendx = dx*σ*Δt*RKs
    tendy = dy*Δt*RKs
    tendz = dz*Δt*RKs

    # pack and return
    return cat(dims=1,x',y',z',dx',zs',ρz',xρz',dy',ys',xy',βz',dz',tendx',tendy',tendz');
end
