include("/home/kloewer/git/Juls.jl/src/Juls.jl")
using .Juls
using BFloat16s
using SoftPosit

Base.round(x::BFloat16, r::RoundingMode{:Up}) = BFloat16(ceil(Float32(x)))
Base.round(x::BFloat16, r::RoundingMode{:Down}) = BFloat16(floor(Float32(x)))
Base.round(x::BFloat16, r::RoundingMode{:Nearest}) = BFloat16(round(Float32(x)))
Base.Int64(x::BFloat16) = Int64(Float32(x))

RunJuls(Float32,
        nx=800,
        output=true,
        Ndays=500.0,
        outpath="/network/aopp/chaos/pred/kloewer/julsdata/forecast/nx800/",
        output_vars=["u","v","η"])

RunJuls(Float64,
        output=true,
        Ndays=50,
        nx=800,
        outpath="/network/aopp/chaos/pred/kloewer/julsdata/forecast/nx800/",
        output_vars=["sst"],
        initial_cond="ncfile",
        init_run_id=0)

RunJuls(Posit16,
        output=true,
        Ndays=50,
        nx=800,
        outpath="/network/aopp/chaos/pred/kloewer/julsdata/forecast/nx800/",
        output_vars=["sst"],
        initial_cond="ncfile",
        init_run_id=0)

RunJuls(Float16,
        output=true,
        Ndays=50,
        nx=800,
        outpath="/network/aopp/chaos/pred/kloewer/julsdata/forecast/nx800/",
        output_vars=["sst"],
        initial_cond="ncfile",
        init_run_id=0)

RunJuls(Float16,
        Tprog=Float32,
        output=true,
        Ndays=50,
        nx=800,
        outpath="/network/aopp/chaos/pred/kloewer/julsdata/forecast/nx800/",
        output_vars=["sst"],
        initial_cond="ncfile",
        init_run_id=0)

RunJuls(BFloat16,
        Tprog=Float32,
        output=true,
        Ndays=50,
        nx=800,
        outpath="/network/aopp/chaos/pred/kloewer/julsdata/forecast/nx800/",
        output_vars=["sst"],
        initial_cond="ncfile",
        init_run_id=0)
