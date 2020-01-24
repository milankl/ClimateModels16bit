include("/home/kloewer/git/ShallowWaters.jl/src/ShallowWaters.jl")
using .ShallowWaters
using BFloat16s
using SoftPosit

Base.round(x::BFloat16, r::RoundingMode{:Up}) = BFloat16(ceil(Float32(x)))
Base.round(x::BFloat16, r::RoundingMode{:Down}) = BFloat16(floor(Float32(x)))
Base.round(x::BFloat16, r::RoundingMode{:Nearest}) = BFloat16(round(Float32(x)))
Base.Int64(x::BFloat16) = Int64(Float32(x))

RunModel(Float32,
        nx=1000,
        L_ratio=4,
        output=true,
        Ndays=100.0,
        outpath="/network/aopp/chaos/pred/kloewer/julsdata/forecast2/nx1000/",
        output_vars=["u","v","η"])
#
# RunModel(Float64,
#         output=true,
#         Ndays=50,
#         nx=400,
#         outpath="/network/aopp/chaos/pred/kloewer/julsdata/forecast2/nx400/",
#         output_vars=["sst"],
#         initial_cond="ncfile",
#         init_run_id=0)
#
# RunModel(Posit16,
#         output=true,
#         Ndays=50,
#         nx=400,
#         outpath="/network/aopp/chaos/pred/kloewer/julsdata/forecast2/nx400/",
#         output_vars=["sst"],
#         initial_cond="ncfile",
#         init_run_id=0)
#
# RunModel(Posit16_2,
#         output=true,
#         Ndays=50,
#         nx=400,
#         outpath="/network/aopp/chaos/pred/kloewer/julsdata/forecast2/nx400/",
#         output_vars=["sst"],
#         initial_cond="ncfile",
#         init_run_id=0)
#
# RunModel(Float16,
#         output=true,
#         Ndays=50,
#         nx=400,
#         outpath="/network/aopp/chaos/pred/kloewer/julsdata/forecast2/nx400/",
#         output_vars=["sst"],
#         initial_cond="ncfile",
#         init_run_id=0)
#
# RunModel(Float16,
#         Tprog=Float32,
#         output=true,
#         Ndays=50,
#         nx=400,
#         outpath="/network/aopp/chaos/pred/kloewer/julsdata/forecast2/nx400/",
#         output_vars=["sst"],
#         initial_cond="ncfile",
#         init_run_id=0)
#
# RunModel(BFloat16,
#         Tprog=Float32,
#         output=true,
#         Ndays=50,
#         nx=400,
#         outpath="/network/aopp/chaos/pred/kloewer/julsdata/forecast2/nx400/",
#         output_vars=["sst"],
#         initial_cond="ncfile",
#         init_run_id=0)
