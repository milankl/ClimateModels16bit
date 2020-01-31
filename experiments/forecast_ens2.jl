include("/home/kloewer/git/ShallowWaters.jl/src/ShallowWaters.jl")
using .ShallowWaters
using FileIO
using BFloat16s
using SoftPosit

Base.round(x::BFloat16, r::RoundingMode{:Up}) = BFloat16(ceil(Float32(x)))
Base.round(x::BFloat16, r::RoundingMode{:Down}) = BFloat16(floor(Float32(x)))
Base.round(x::BFloat16, r::RoundingMode{:Nearest}) = BFloat16(round(Float32(x)))
Base.Int64(x::BFloat16) = Int64(Float32(x))


"""Finds the first gap in a list of integers."""
function gap(a::Array{Int,1})
    try
        return minimum([i for i in minimum(a):maximum(a) if ~(i in a)])
    catch
        return maximum(a)+1
    end
end

function get_run_id(path::String,order::String="continue")
    runlist = filter(x->startswith(x,"run"),readdir(path))
    existing_runs = [parse(Int,id[4:end]) for id in runlist]
    if length(existing_runs) == 0           # if no runfolder exists yet
        run_id = 0
    else                                    # create next folder
        if order == "fill"  # find the smallest gap in runfolders
            run_id = gap(existing_runs)
        elseif order == "continue" # find largest folder and count one up
            run_id = maximum(existing_runs)+1
        end
    end
    return run_id
end

path = "/network/aopp/chaos/pred/kloewer/julsdata/forecast2/"
startis = load(joinpath(path,"starti.jld2"))["starti"]
outpath = joinpath(path,"comm","Float64")

for i in 1:100
    run_id = get_run_id(outpath,"fill")

    starti = startis[run_id+1]

    RunModel(Float64,
            output=true,
            Ndays=300.0,
            output_dt=12,
            outpath=outpath,
            initial_cond="ncfile",
            output_vars=["u","η","sst"],
            initpath=path,
            init_run_id=1,
            init_starti=starti,
            get_id_mode="fill")

end
