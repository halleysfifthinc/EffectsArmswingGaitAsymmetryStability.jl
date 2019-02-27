module SteadyStateStability

using Biomechanics, MAT, ChaosTools, Interpolations, UnsafeArrays, DSP

# Standard library
using Statistics, DelimitedFiles, LinearAlgebra

export readsegment,
       readsstrials,
       analyzetrial

export DSSteadyState,
       SteadyStateSeg

abstract type DSSteadyState <: AbstractDataSource end

struct SteadyStateSeg <: DSSteadyState
    events::Dict{Symbol,Vector{Float64}}
    data::Dict{Symbol,Matrix{Float64}}
end

SteadyStateSeg() = SteadyStateSeg(Dict{Symbol,Matrix}(), Dict{Symbol,AbstractVector}())

function readsegment(DSData::Type{<:AbstractDataSource},
                     trial::Trial{<:AbstractDataSource},
                     st::Float64 = 0.0,
                     en::Float64 = Inf;
                     events::Vector = [],
                     ts::Vector = [],
                     fs = 100)
    isfile(trial.path) || throw(ArgumentError("trial $(trial.path) exist"))
    st >= 0.0 || throw(ArgumentError("start time must be positive"))
    st <= en || throw(ArgumentError("end time must be greater than start time"))

    evnames = Symbol.(events)
    file = matopen(trial.path)

    revents = Dict{Symbol,Vector{Float64}}()
    for e in events
        if exists(file, string(e))
            syme = Symbol(e)
            tmp = read(file, string(e))[1]
            if tmp isa AbstractArray
                revents[syme] = vec(tmp)
            else
                revents[syme] = [tmp]
            end
            strt = findfirst(x -> x >= st, revents[syme])
            if strt == 0
                # @warn "no $e events during given start and end times"
                delete!(revents, syme)
                break
            end
            endi = findlast(x -> x <= en, revents[syme])
            revents[syme] = revents[syme][strt:endi] .- st .+ (1/fs) # Shift events to be index accurate for data subsection
        else
            # @warn "Requested event $e does not exist in source data"
        end
    end

    data = Dict{Symbol,Matrix{Float64}}()
    for t in ts
        if exists(file, string(t))
            symt = Symbol(t)
            data[symt] = read(file, string(t))[1]
            len = size(data[symt], 1)
            strti = round(Int, st*fs)
            endi = en == Inf ? len : min(len, round(Int, en*fs))
            data[symt] = data[symt][strti:endi, :]
        else
            @warn "Requested time series $t does not exist in source data"
        end
    end

    return Segment(trial, Dict{Symbol,Any}(), DSData(revents, data))
end

const fs = 100

function readsstrials(rootdir::String, subs=1:15)
    genV3D = joinpath(rootdir, "data", "generated", "ARMS_STAB")

    trialnames = [  sym*"-"*arms for sym in ["sym", "asym"] for arms in ["none", "norm", "excess"] ].*".mat"

    sstrials = Vector{Trial}()

    for sub in subs
        subj = lpad(sub, 2, '0')
        pref = joinpath(genV3D, "Subject $subj", "export", "steady-state")
        subtrials = readdir(pref)
        filter!(subtrials) do file
            basename(file) ∈ trialnames
        end

        subtrials .= joinpath.(pref, subtrials)

        for trial in subtrials
            conds = Dict{Symbol, Symbol}()
            m = match(r"[\\,\/](?<sym>(sym|asym))-(?<arms>(none|norm|excess))", trial)
            conds[:sym] = Symbol(m[:sym])
            conds[:arms] = Symbol(m[:arms])

            push!(sstrials, Trial{DSSteadyState}(sub, splitext(basename(trial))[1], trial, conds))
        end
    end

    return sstrials
end

# Lamb and Stöckl 2014 doi:10.1016/j.clinbiomech.2014.03.008
function continuousphase(θ, events::AbstractVector{<:Integer})
    mi, ma = avgextrema(θ, events)
    θcent = θ .- Ref(mi) .- Ref((ma - mi)/2)
    Hθ = hilbert(θcent)

    return angle.(Hθ)
end

# Minor type-piracy
Base.:*(::Nothing, ::Number) = nothing
Base.:-(::Nothing, ::Number) = nothing
Base.min(::Nothing, x::Number) = x
Base.min(x::Number, ::Nothing) = x

function rom(x::AbstractVector{T}, events::AbstractVector{<:Integer}) where {T}
    l = length(events) - 1
    E = Array{T}(undef, l, 2)
    _rangedextrema!(E, x, events, l)

    roms = diff(E; dims=2)

    avgrom = mean(roms)
    stdrom = std(roms)

    (avgrom, stdrom)
end

function avgextrema(x::AbstractVector{T}, events::AbstractVector{<:Integer}) where {T}
    l = length(events) - 1
    E = Array{T}(undef, l, 2)
    _rangedextrema!(E, x, events, l)

    # ma = mean(@view(E[:,2]))
    # mi = mean(@view(E[:,1]))

    avgex = mean(E; dims=1)

    (avgex[1], avgex[2])
end

function _rangedextrema!(E, x, events, l)
    @inbounds for i in 1:l
        GC.@preserve x begin
            E[i,1], E[i,2] = extrema(uview(x, events[i]:(events[i+1]-1)))
        end
    end
    return E
end

# Asymmetry function (see Plotnik et al. 2005 doi:10.1002/ana.20452)
asymmetry(l,r) = abs(log(min(l,r)/max(l,r)))*100

function analyzetrial(trial::Trial, numstrides::Int)
    cols = [ "TrunkLinVel", "TrunkAngVel", "RFootPos", "LFootPos", "LHip", "RHip", "LShoulder", "RShoulder" ]
    seg = readsegment(SteadyStateSeg, trial, 25.0; events=[:RFC, :LFC], ts=cols)
    insuffstrerr = "insufficient number of strides in trial $(trial.name), subject $(trial.subject)"
    length(seg.data.events[:RFC]) < numstrides+1 && throw(ArgumentError(insuffstrerr))

    results = Dict{Symbol,Any}()

    ########################################
    ## Interleaving steps
    ########################################

    # Gather rfc and lfc info into a NamedTuple, with named fields `time`, `pos`, and `leg`, (`:RFC` or `:LFC`)
    rfc = [ (time=seg.data.events[:RFC][i],
             pos=SVector{2}(seg.data.data[:RFootPos][
                    round(Int, seg.data.events[:RFC][i]*fs)
                    ,1:2]),
             leg=:RFC) for i in 1:(numstrides+1) ]
    if seg.data.events[:LFC][1] < seg.data.events[:RFC][1]
        lfc = [ (time=seg.data.events[:LFC][i],
                 pos=SVector{2}(seg.data.data[:LFootPos][
                        round(Int, seg.data.events[:LFC][i]*fs)
                        ,1:2]),
                 leg=:LFC) for i in 2:(numstrides+2) ]
    else
        lfc = [ (time=seg.data.events[:LFC][i],
                 pos=SVector{2}(seg.data.data[:LFootPos][
                        round(Int, seg.data.events[:LFC][i]*fs)
                        ,1:2]),
                 leg=:LFC) for i in 1:(numstrides+1) ]
    end

    # Sort all steps and check that they are alternating--no double steps--and that they are in the expected locations:
    # Odd elements need to be RFC for later assumptions made on whether the odd/even elements of the diff are right or left steps
    steps = sort([rfc; lfc], by=(x -> x.time))

    while any(x -> x.leg === :LFC, steps[isodd.(axes(steps, 1))])
        badstep = min(findfirst(x -> x.leg === :LFC, steps[isodd.(axes(steps, 1))])*2-1,
                      findfirst(x -> x.leg === :RFC, steps[iseven.(axes(steps, 1))])*2)
        steps = deleteat!(steps, badstep)
    end
    all(x -> x.leg === :RFC, steps[isodd.(axes(steps, 1))]) || throw(DomainError("Odd elements aren't all RFC's"))
    all(x -> x.leg === :LFC, steps[iseven.(axes(steps, 1))]) || throw(DomainError("Even elements aren't all LFC's"))

    ########################################
    ## Gait variability
    ########################################

    # After `diff`ing, odd elements will be left steps, even elements will be right steps

    steptimes = diff(getindex.(steps, :time))
    if trial.conds[:sym] == :asym
        # The speed of the treadmill belt for the planted foot is added to the distance between successive steps
        # e.g. For the left step, the right foot is planted, therefore the effective AP distance traveled by the
        # left foot is speed of the right treadmill belt plus the absolute distance between the right and left
        # footstrikes
        stepcoords = diff(getindex.(steps, :pos)) .+
                        ( isodd(i) ? SVector{2}(0.0, 0.96*steptimes[i]) : SVector{2}(0.0, 1.2*steptimes[i])
                         for i in 1:(length(steps)-1) )
    else
        stepcoords = diff(getindex.(steps, :pos)) .+
                        ( SVector{2}(0.0, 1.2*steptimes[i]) for i in 1:(length(steps)-1) )
    end

    # Spatial and temporal step descriptives
    left_stepcoords = stepcoords[isodd.(axes(stepcoords,1))]
    left_steplavg = mean(abs.(getindex.(left_stepcoords,2)))
    left_steplstd = std(abs.(getindex.(left_stepcoords,2)))

    right_stepcoords = stepcoords[iseven.(axes(stepcoords,1))]
    right_steplavg = mean(abs.(getindex.(right_stepcoords,2)))
    right_steplstd = std(abs.(getindex.(right_stepcoords,2)))

    # Push into results Dict
    results[:left_steplavg] = left_steplavg
    results[:left_steplstd] = left_steplstd
    results[:right_steplavg] = right_steplavg
    results[:right_steplstd] = right_steplstd

    left_avgsteptime = mean(steptimes[isodd.(axes(steptimes,1))])
    right_avgsteptime = mean(steptimes[iseven.(axes(steptimes,1))])

    # Asymmetry indices
    results[:stepasym_spatial] = asymmetry(left_steplavg, right_steplavg)
    results[:stepasym_temporal] = asymmetry(left_avgsteptime, right_avgsteptime)

    # Step width isn't confounded by induced asymmetric gait
    results[:stepwavg] = mean(abs.(getindex.(stepcoords, 1)))
    results[:stepwstd] = std(abs.(getindex.(stepcoords, 1)))

    ########################################
    # Arm swing
    ########################################
    rfc = round.(Int, seg.data.events[:RFC][1:numstrides+1] .* fs) .- 1 # Subtract 1 because the θ and ω have the first element trimmed

    # Drop last element (it is a NaN)
    rhip_θ = seg.data.data[:RHip][1:end-1,1]
    lhip_θ = seg.data.data[:LHip][1:end-1,1]
    rsho_θ = seg.data.data[:RShoulder][1:end-1,1]
    lsho_θ = seg.data.data[:LShoulder][1:end-1,1]

    # ROM and swing asymmetry
    rsho_avgrom, rsho_stdrom = rom(rsho_θ, rfc)
    lsho_avgrom, lsho_stdrom = rom(lsho_θ, rfc)
    results[:swingasym] = asymmetry(lsho_avgrom, rsho_avgrom)

    results[:rsho_avgrom] = rsho_avgrom
    results[:lsho_avgrom] = lsho_avgrom
    results[:rsho_stdrom] = rsho_stdrom
    results[:lsho_stdrom] = lsho_stdrom

    # Continuous phase
    rhip_crp = continuousphase(rhip_θ, rfc)
    lhip_crp = continuousphase(lhip_θ, rfc)
    rsho_crp = continuousphase(rsho_θ, rfc)
    lsho_crp = continuousphase(lsho_θ, rfc)

    # The subtraction of the two angles has the possibility of introducing discontinuties (which
    # will interfere with the interpolation of the time normalizing)
    crp_lsho_rhip = Unwrap.unwrap(rhip_crp; range=2pi) - Unwrap.unwrap(lsho_crp; range=2pi)
    crp_rsho_lhip = Unwrap.unwrap(lhip_crp; range=2pi) - Unwrap.unwrap(rsho_crp; range=2pi)

    norm_crp_lsho_rhip = timenormalize(crp_lsho_rhip, rfc)
    norm_crp_rsho_lhip = timenormalize(crp_rsho_lhip, rfc)

    # Circular SD (Fisher 1993)
    S = sum(sin.(reshape(norm_crp_lsho_rhip, 100, numstrides)); dims=2)./numstrides
    C = sum(cos.(reshape(norm_crp_lsho_rhip, 100, numstrides)); dims=2)./numstrides
    sdcrp_lsho_rhip = map((c, s) -> sqrt(-2*log(sqrt(c^2 + s^2))), C, S)

    # Circular mean (Fisher 1993)
    s = sum(sin.(sdcrp_lsho_rhip))/numstrides
    c = sum(cos.(sdcrp_lsho_rhip))/numstrides
    msdcrp_lsho_rhip = atan(s, c)

    S = sum(sin.(reshape(norm_crp_rsho_lhip, 100, numstrides)); dims=2)./numstrides
    C = sum(cos.(reshape(norm_crp_rsho_lhip, 100, numstrides)); dims=2)./numstrides
    sdcrp_rsho_lhip = map((c, s) -> sqrt(-2*log(sqrt(c^2 + s^2))), C, S)

    s = sum(sin.(sdcrp_rsho_lhip ))/numstrides
    c = sum(cos.(sdcrp_rsho_lhip ))/numstrides
    msdcrp_rsho_lhip = atan(s, c)

    results[:msdcrp_lsho_rhip] = rad2deg(msdcrp_lsho_rhip)
    results[:msdcrp_rsho_lhip] = rad2deg(msdcrp_rsho_lhip)

    ########################################
    ## Lyapunov exponent
    ########################################
    frames = round(Int, seg.data.events[:RFC][1]*fs):round(Int, seg.data.events[:RFC][numstrides+1]*fs)

    data = Matrix{Float64}(undef, length(frames), 2*3)
    three = 1:3
    for (i, t) in enumerate(["TrunkLinVel", "TrunkAngVel"])
        data[:,three .+ (i-1)*3] = seg.data.data[Symbol(t)][frames,:]
    end

    NW = Matrix{Float64}(undef, 12_500, 2*3)
    for i in 1:size(NW,2)
        itp = interpolate(data[:,i], BSpline(Cubic(Line(Interpolations.OnGrid()))))
        NW[:,i] .= itp(range(1, stop=length(frames), length=12_500))
    end

    linvelstd = sum(std(view(NW, :, 1:3); dims=1))
    angvelstd = sum(std(view(NW, :, 4:6); dims=1))

    NW[:,1:3] ./= linvelstd
    NW[:,4:6] ./= angvelstd

    d = Dataset(NW)
    r = reconstruct(d, 2, 25)

    ks = 0:50
    E = numericallyapunov(r, ks; w=100, ntype=FixedMassNeighborhood(1), distance=Euclidean())

    λ = ([ ones(length(ks)) ks./100 ] \ E)[2]
    results[:lambdaS] = λ

    return AnalyzedSegment(seg, results)
end

end
