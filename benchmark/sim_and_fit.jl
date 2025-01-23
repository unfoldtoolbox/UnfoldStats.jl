using UnfoldSim
using Unfold, UnfoldMixedModels
using StableRNGs
using DataFrames

function define_simulation(sim_type, β, σs; n_subjects = 30, n_items = 100, noiselevel = 2)

    # Create design
    conditions = get_conditions(sim_type)

    design = MultiSubjectDesign(;
        n_subjects = n_subjects,
        n_items = n_items,
        items_between = conditions,
    )

    # Specify component
    basis = p100()
    formula = get_formula(sim_type)
    signal = create_component(sim_type, basis, formula, β, σs)

    # Specify inter-onset distribution
    onset = UniformOnset(; width = 50, offset = 1)

    # Specify noise
    noise = PinkNoise(; noiselevel = noiselevel)

    return Simulation(design, signal, onset, noise)
end

function sim_and_fit(
    sim_type::SimulationType,
    simulation::Simulation,
    model_type::Type{<:UnfoldModel};
    seed::Integer = 1,
)

    # At the moment, the function is not implemented for mixed models
    #ext = Base.get_extension(Unfold, :UnfoldMixedModelsExt)
    if model_type in [UnfoldLinearMixedModel, UnfoldLinearMixedModelContinuousTime]
        throw("Not implemented.")
    end

    # Set parameter(s) for data simulation
    if model_type == UnfoldLinearModel
        return_epoched = true
    else # UnfoldLinearModelContinuousTime
        return_epoched = false
    end

    # Simulate data
    data, events = simulate_data(sim_type, simulation, return_epoched, seed)

    # Create event dict containing basis function(s)/times and formula(s) for all events
    event_dict = create_event_dict(sim_type, model_type, simulation)

    # Fit an Unfold model for each subject
    subject_list = unique(events.subject)
    model_list = model_type[]

    # Slice the data by its last dimension (i.e. the subject dimension)
    data_slices = eachslice(data, dims = ndims(data))

    for s = 1:size(data, ndims(data))
        m = fit(
            UnfoldModel,
            event_dict,
            subset(events, :subject => ByRow(==(subject_list[s]))),
            data_slices[s],
        )
        push!(model_list, m)
    end

    models = DataFrame(subject = subject_list, unfoldmodel = model_list)

    return models
end

simulate_data(sim_type::T, simulation, return_epoched, seed) where {T} =
    simulate_data(EventStyle(T), simulation, return_epoched, seed)

function simulate_data(::SingleEventType, simulation, return_epoched, seed)

    # Simulate data
    data, events = simulate(StableRNG(seed), simulation; return_epoched = return_epoched)

    return data, events
end

function simulate_data(::MultipleEventTypes, simulation, return_epoched, seed)

    # Simulate data
    data, events = simulate(StableRNG(seed), simulation; return_epoched = return_epoched)

    # Add an event column to the events df and assign each event to an event type
    events[!, :event] = repeat(["stim", "fix"], size(events, 1) ÷ 2)

    return data, events
end

create_event_dict(sim_type::T, model_type, simulation) where {T} = create_event_dict(
    EventStyle(T),
    PredictorStyle(T),
    model_type::Type{<:UnfoldModel},
    simulation,
)

function create_event_dict(
    ::MultipleEventTypes,
    ::ManyPredictors,
    model_type::Type{<:UnfoldModel},
    simulation,
)
    # Create times vector/basis function(s) (for model fitting)
    if model_type == UnfoldLinearModel
        #times = axes(data, 1)
        times = 1:UnfoldSim.maxlength(simulation.components)
        t_stim = times
        t_fix = times
    else # UnfoldLinearModelContinuousTime
        t_stim = firbasis(τ = (-0.1, 1.5), sfreq = 100, name = "stim")
        t_fix = firbasis(τ = (-0.1, 1), sfreq = 100, name = "fix")
    end

    # Define formula(s)
    f_stim = @formula 0 ~ 1 + continuous
    f_fix = @formula 0 ~ 1 + spl(continuous, 4) + continuous + condition * pet

    # Combine basis function(s)/times and formula(s) with the corresponding event
    event_dict = Dict("stim" => (f_stim, t_stim), "fix" => (f_fix, t_fix))

    return event_dict
end

function create_event_dict(
    ::MultipleEventTypes,
    ::OnlySplines,
    model_type::Type{<:UnfoldModel},
    simulation,
)
    # Create times vector/basis function(s) (for model fitting)
    if model_type == UnfoldLinearModel
        #times = axes(data, 1)
        times = 1:UnfoldSim.maxlength(simulation.components)
        t_stim = times
        t_fix = times
    else # UnfoldLinearModelContinuousTime
        t_stim = firbasis(τ = (-0.1, 1.5), sfreq = 100, name = "stim")
        t_fix = firbasis(τ = (-0.1, 1), sfreq = 100, name = "fix")
    end

    # Define formula(s)
    f_stim = @formula 0 ~ 1
    f_fix = @formula 0 ~ 1 + spl(continuous, 4)

    # Combine basis function(s)/times and formula(s) with the corresponding event
    event_vec = ["stim" => (f_stim, t_stim), "fix" => (f_fix, t_fix)]

    return event_vec
end

function create_event_dict(
    ::SingleEventType,
    ::ManyPredictors,
    model_type::Type{<:UnfoldModel},
    simulation,
)
    # Create times vector/basis function(s) (for model fitting)
    if model_type == UnfoldLinearModel
        #times = axes(data, 1)
        t = 1:UnfoldSim.maxlength(simulation.components)
    else # UnfoldLinearModelContinuousTime
        t = firbasis((-0.1, 1.0), 100)
    end

    # Define formula(s)
    f = @formula 0 ~ 1 + spl(continuous, 4) + continuous + condition * pet

    # Combine basis function(s)/times and formula(s) with the corresponding event
    event_vec = [Any => (f, t)]

    return event_vec
end

function create_event_dict(
    ::SingleEventType,
    ::OnlySplines,
    model_type::Type{<:UnfoldModel},
    simulation,
)
    # Create times vector/basis function(s) (for model fitting)
    if model_type == UnfoldLinearModel
        #times = axes(data, 1)
        t = 1:UnfoldSim.maxlength(simulation.components)
    else # UnfoldLinearModelContinuousTime
        t = firbasis((-0.1, 1.0), 100)
    end

    # Define formula(s)
    f = @formula 0 ~ 1 + spl(continuous, 4)

    # Combine basis function(s)/times and formula(s) with the corresponding event
    event_vec = [Any => (f, t)]

    return event_vec
end

get_conditions(sim_type::T) where {T} = get_conditions(PredictorStyle(T))

function get_conditions(::OnlySplines)
    conditions = Dict(:continuous => range(-5, 5, length = 50))
    return conditions
end

function get_conditions(::ManyPredictors)
    conditions = Dict(
        :continuous => range(-5, 5, length = 25),
        :condition => ["face", "car"],
        :pet => ["cat", "dog"],
    )
    return conditions
end

get_formula(sim_type::T) where {T} = get_formula(PredictorStyle(T))

function get_formula(::OnlySplines)
    formula = @formula(
        0 ~
            1 +
            continuous +
            continuous^2 +
            continuous^3 +
            (1 + continuous + continuous^2 + continuous^3 | subject)
    )
    return formula
end

function get_formula(::ManyPredictors)
    formula = @formula(
        0 ~
            1 +
            continuous +
            continuous^2 +
            continuous^3 +
            condition +
            pet +
            (1 + continuous + continuous^2 + continuous^3 | subject)
    )
    return formula
end

create_component(sim_type::T, basis, formula, β, σs) where {T} =
    create_component(ChannelStyle(T), basis, formula, β, σs)

function create_component(::SingleChannel, basis, formula, β, σs)
    signal = MixedModelComponent(; basis = basis, formula = formula, β = β, σs = σs)
    return signal
end

function create_component(::MultiChannel, basis, formula, β, σs)
    signal = MixedModelComponent(; basis = basis, formula = formula, β = β, σs = σs)
    # Wrap the component in a multichannel component
    # Load headmodel
    hart = headmodel(type = "hartmut")
    source_idx = UnfoldSim.closest_src(hart, "Left Postcentral Gyrus")
    projection = UnfoldSim.magnitude(hart)
    # Only use the first channels/electrodes
    channels = 1:3

    multichannel_signal =
        MultichannelComponent(signal, projection[channels, source_idx], NoNoise())
    return multichannel_signal
end
