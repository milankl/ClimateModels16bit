include("/home/kloewer/git/ShallowWaters.jl/src/ShallowWaters.jl")
using .ShallowWaters

RunModel(Float32,
        output=true,
        Ndays=1000.0,
        outpath="/network/aopp/chaos/pred/kloewer/julsdata/forecast",
        output_vars=["u","v","η"])

RunModel(Float32,
        output=true,
        Ndays=50*365,
        outpath="/network/aopp/chaos/pred/kloewer/julsdata/forecast",
        output_vars=["u","v","η"],
        initial_cond="ncfile",
        init_run_id=0)
