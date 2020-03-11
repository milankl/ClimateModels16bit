include("/home/kloewer/git/ShallowWaters.jl/src/ShallowWaters.jl")
using .ShallowWaters
using BFloat16s
using SoftPosit
#
# Base.round(x::BFloat16, r::RoundingMode{:Up}) = BFloat16(ceil(Float32(x)))
# Base.round(x::BFloat16, r::RoundingMode{:Down}) = BFloat16(floor(Float32(x)))
# Base.round(x::BFloat16, r::RoundingMode{:Nearest}) = BFloat16(round(Float32(x)))
# Base.Int64(x::BFloat16) = Int64(Float32(x))

# RunModel(Float32,
#         nx=400,
#         output=true,
#         Ndays=500.0,
#         outpath="/network/aopp/chaos/pred/kloewer/julsdata/forecast2/nx400comm/",
#         output_vars=["u","v","η"])

# RunModel(Float64,
#         output=true,
#         Ndays=50,
#         nx=400,
#         outpath="/network/aopp/chaos/pred/kloewer/julsdata/forecast2/nx400comm/",
#         output_vars=["sst"],
#         initial_cond="ncfile",
#         init_run_id=0)

# RunModel(Float64,
#         Tcomm=Float16,
#         output=true,
#         Ndays=50,
#         nx=400,
#         outpath="/network/aopp/chaos/pred/kloewer/julsdata/forecast2/nx400comm/",
#         output_vars=["sst"],
#         initial_cond="ncfile",
#         init_run_id=0)
#
# RunModel(Float64,
#         Tcomm=BFloat16,
#         output=true,
#         Ndays=50,
#         nx=400,
#         outpath="/network/aopp/chaos/pred/kloewer/julsdata/forecast2/nx400comm/",
#         output_vars=["sst"],
#         initial_cond="ncfile",
#         init_run_id=0)
#
# RunModel(Float64,
#         Tcomm=Posit16,
#         output=true,
#         Ndays=50,
#         nx=400,
#         outpath="/network/aopp/chaos/pred/kloewer/julsdata/forecast2/nx400comm/",
#         output_vars=["sst"],
#         initial_cond="ncfile",
#         init_run_id=0)
#
# RunModel(Float64,
#         Tcomm=Posit16_2,
#         output=true,
#         Ndays=50,
#         nx=400,
#         outpath="/network/aopp/chaos/pred/kloewer/julsdata/forecast2/nx400comm/",
#         output_vars=["sst"],
#         initial_cond="ncfile",
#         init_run_id=0)

RunModel(BFloat16sr,
        output=true,
        Ndays=50,
        nx=400,
        outpath="/network/aopp/chaos/pred/kloewer/julsdata/forecast/nx400/",
        output_vars=["sst"],
        initial_cond="ncfile",
        init_run_id=0)
