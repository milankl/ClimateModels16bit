include("/home/kloewer/git/ShallowWaters.jl/src/ShallowWaters.jl")
using .ShallowWaters
using BFloat16s

Base.round(x::BFloat16, r::RoundingMode{:Up}) = BFloat16(ceil(Float32(x)))
Base.round(x::BFloat16, r::RoundingMode{:Down}) = BFloat16(floor(Float32(x)))
Base.round(x::BFloat16, r::RoundingMode{:Nearest}) = BFloat16(round(Float32(x)))
Base.Int64(x::BFloat16) = Int64(Float32(x))

RunModel(Float32,
        nx=2000,
        L_ratio=4,
        cfl=0.9,
        output=true,
        Ndays=500.0,
        g=1e-1,
        output_dt=48,
        outpath="/network/aopp/chaos/pred/kloewer/julsdata/forecast2/nx2000/",
        output_vars=["u","v","η"])

RunModel(BFloat16,
        Tprog=Float32,
        L_ratio=4,
        output=true,
        Ndays=100,
        cfl=0.9,
        nx=2000,
        g=1e-1,
        outpath="/network/aopp/chaos/pred/kloewer/julsdata/forecast2/nx2000/",
        output_vars=["sst"],
        initial_cond="ncfile",
        init_run_id=0)
