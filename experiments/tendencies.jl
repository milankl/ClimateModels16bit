include("/home/kloewer/git/ShallowWaters.jl/src/ShallowWaters.jl")
using .ShallowWaters
using SoftPosit
using BFloat16s

Base.round(x::BFloat16, r::RoundingMode{:Up}) = BFloat16(ceil(Float32(x)))
Base.round(x::BFloat16, r::RoundingMode{:Down}) = BFloat16(floor(Float32(x)))
Base.round(x::BFloat16, r::RoundingMode{:Nearest}) = BFloat16(round(Float32(x)))
Base.Int64(x::BFloat16) = Int64(Float32(x))

# RunModel(Float64,
#         output=true,
#         Ndays=50,
#         nx=400,
#         initpath="/network/aopp/chaos/pred/kloewer/julsdata/forecast/nx400/",
#         outpath="/network/aopp/chaos/pred/kloewer/julsdata/forecast/tend/",
#         output_vars=["u","v","η","du","dv","dη"],
#         initial_cond="ncfile",
#         init_run_id=0)
#
# RunModel(Float16,
#         output=true,
#         Ndays=50,
#         nx=400,
#         initpath="/network/aopp/chaos/pred/kloewer/julsdata/forecast/nx400/",
#         outpath="/network/aopp/chaos/pred/kloewer/julsdata/forecast/tend/",
#         output_vars=["u","v","η","du","dv","dη"],
#         initial_cond="ncfile",
#         init_run_id=0)
#
# RunModel(Float16,
#         Tprog=Float32,
#         output=true,
#         Ndays=50,
#         nx=400,
#         initpath="/network/aopp/chaos/pred/kloewer/julsdata/forecast/nx400/",
#         outpath="/network/aopp/chaos/pred/kloewer/julsdata/forecast/tend/",
#         output_vars=["u","v","η","du","dv","dη"],
#         initial_cond="ncfile",
#         init_run_id=0)

# RunModel(BFloat16,
#         Tprog=Float32,
#         output=true,
#         Ndays=50,
#         nx=400,
#         initpath="/network/aopp/chaos/pred/kloewer/julsdata/forecast/nx400/",
#         outpath="/network/aopp/chaos/pred/kloewer/julsdata/forecast/tend/",
#         output_vars=["u","v","η","du","dv","dη"],
#         initial_cond="ncfile",
#         init_run_id=0)

# RunModel(Posit16,
#         output=true,
#         Ndays=50,
#         nx=400,
#         initpath="/network/aopp/chaos/pred/kloewer/julsdata/forecast/nx400/",
#         outpath="/network/aopp/chaos/pred/kloewer/julsdata/forecast/tend/",
#         output_vars=["u","v","η","du","dv","dη"],
#         initial_cond="ncfile",
#         init_run_id=0)
#
RunModel(Posit16_2,
        output=true,
        Ndays=50,
        nx=400,
        initpath="/network/aopp/chaos/pred/kloewer/julsdata/forecast/nx400/",
        outpath="/network/aopp/chaos/pred/kloewer/julsdata/forecast/tend/",
        output_vars=["u","v","η","du","dv","dη"],
        initial_cond="ncfile",
        init_run_id=0)
