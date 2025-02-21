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
    formula = get_sim_formula(sim_type)
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
    if model_type in [UnfoldLinearMixedModelContinuousTime]
        throw("Not implemented.")
    end

    # Set parameter(s) for data simulation
    if model_type == UnfoldLinearModel || model_type == UnfoldLinearMixedModel
        return_epoched = true
    else # UnfoldLinearModelContinuousTime
        return_epoched = false
    end

    # Simulate data
    data, events = simulate_data(sim_type, simulation, return_epoched, seed)

    # Create event dict containing basis function(s)/times and formula(s) for all events
    event_vec = create_event_vector(sim_type, model_type, simulation)

    if model_type == UnfoldLinearMixedModel
        if length(size(data)) == 4
            # if channels exist already
            data = reshape(data, size(data)[1:end-2]..., :)
        else
            data = reshape(data, 1, size(data)[1:end-2]..., :)
        end
        m = fit(UnfoldModel, event_vec, events, data)
        return DataFrame(subject = ["LMM - all subjects"], unfoldmodel = [m], data = [data])

    else
        # Fit an Unfold model for each subject
        subject_list = unique(events.subject)
        model_list = model_type[]

        # Slice the data by its last dimension (i.e. the subject dimension)
        data_slices = eachslice(data, dims = ndims(data))

        for s = 1:size(data, ndims(data))
            m = fit(
                UnfoldModel,
                event_vec,
                subset(events, :subject => ByRow(==(subject_list[s]))),
                data_slices[s],
            )
            push!(model_list, m)
        end

        models =
            DataFrame(subject = subject_list, unfoldmodel = model_list, data = data_slices)

        return models
    end
end

simulate_data(sim_type::T, simulation, return_epoched, seed) where {T} =
    simulate_data(EventStyle(T), simulation, return_epoched, seed)

function simulate_data(::SingleEventType, simulation, return_epoched, seed)

    # Simulate data
    data, events =
        UnfoldSim.simulate(StableRNG(seed), simulation; return_epoched = return_epoched)

    return data, events
end

function simulate_data(::MultipleEventTypes, simulation, return_epoched, seed)

    # Simulate data
    data, events =
        UnfoldSim.simulate(StableRNG(seed), simulation; return_epoched = return_epoched)

    # Add an event column to the events df and assign each event to an event type
    events[!, :event] = repeat(["stim", "fix"], size(events, 1) ÷ 2)

    return data, events
end

create_event_vector(sim_type::T, model_type, simulation) where {T} = create_event_vector(
    EventStyle(T),
    PredictorStyle(T),
    model_type::Type{<:UnfoldModel},
    simulation,
)

#=
function create_event_vector(
    ::MultipleEventTypes,
    ::ManyPredictors,
    model_type::Type{<:UnfoldModel},
    simulation,
)
    # Create times vector/basis function(s) (for model fitting)
    t_stim, t_fix = create_timevec(model_type, UnfoldSim.maxlength(simulation.components))

    # Combine basis function(s)/times and formula(s) with the corresponding event
    event_vec = Dict("stim" => (f_stim, t_stim), "fix" => (f_fix, t_fix))

    return event_vec
end
=#
function create_timevec(model_type, event_type, maxlength)
    # Create times vector/basis function(s) (for model fitting)
    if model_type == UnfoldLinearModel || model_type == UnfoldLinearMixedModel
        #times = axes(data, 1)
        times = 1:maxlength
        t_stim = times
        t_fix = times
    else # UnfoldLinearModelContinuousTime
        t_stim = firbasis(τ = (-0.1, 1.5), sfreq = 100)
        t_fix = firbasis(τ = (-0.1, 1), sfreq = 100)
    end
    return t_stim, t_fix

end

create_timevec(model_type, event_type::SingleEventType, ml) =
    create_timevec(model_type, MultipleEventTypes(), ml)[2]

function create_event_vector(
    event_style::EventStyle,
    predictor_style::PredictorStyle,
    model_type::Type{<:UnfoldModel},
    simulation,
)
    times =
        create_timevec(model_type, event_style, UnfoldSim.maxlength(simulation.components))
    # Define formula(s)
    forms = get_fit_formulas(predictor_style, event_style, model_type)

    @debug typeof(times) typeof(forms)

    # Combine basis function(s)/times and formula(s) with the corresponding event
    return create_eventvec(times, forms)

end

create_eventvec(times::Tuple, forms::Vector) =
    ["stim" => (forms[1], times[1]), "fix" => (forms[2], times[2])]
create_eventvec(times, forms) = [Any => (forms, times)]
#get_fit_formulas(predictor_style::ManyPredictors,event_style,model_type) =
#    @formula 0 ~ 1 + spl(continuous, 4) + continuous + condition * pet
#get_fit_formulas(predictor_style::OnlySplines,event_style::SingleEventType,model_type) = @formula 0 ~ 1 + spl(continuous, 4)

get_fit_formulas(
    predictor_style::OnlySplines,
    event_style::MultipleEventTypes,
    model_type,
) = [@formula(0 ~ 1), @formula(0 ~ 1 + spl(continuous, 4))]
get_fit_formulas(
    predictor_style::ManyPredictors,
    event_style::MultipleEventTypes,
    model_type,
) = [
    @formula(0 ~ 1 + continuous),
    @formula(0 ~ 1 + spl(continuous, 4) + continuous + condition * pet),
]
get_fit_formulas(
    predictor_style::CategoricalPredictors,
    event_style::MultipleEventTypes,
    model_type::Type{UnfoldLinearMixedModel},
) = [@formula(0 ~ 1 + continuous + (1 | item)), @formula(0 ~ 1 + pet + (1 + pet | subject))]

get_fit_formulas(predictor_style, event_style::SingleEventType, model_type) =
    get_fit_formulas(predictor_style, MultipleEventTypes(), model_type)[2]




#=
function create_event_vector(
    ::SingleEventType,
    predictorstyle::Union{ManyPredictors,OnlySplines},
    model_type::Type{<:UnfoldModel},
    simulation,
)
    # Create times vector/basis function(s) (for model fitting)
    _, t = create_timevec(model_type, UnfoldSim.maxlength(simulation.components))


    # Define formula(s)
    f = get_fit_formulas(predictor_style)


    # Combine basis function(s)/times and formula(s) with the corresponding event
    event_vec = [Any => (f, t)]

    return event_vec
end

=#
get_conditions(sim_type::T) where {T} = get_conditions(PredictorStyle(T))
get_conditions(::CategoricalPredictors) =
    Dict(:pet => ["cat", "dog"], :condition => ["face", "car"])
get_conditions(::OnlySplines) = Dict(:continuous => range(-5, 5, length = 50))
get_conditions(::ManyPredictors) = Dict(
    :continuous => range(-5, 5, length = 25),
    :condition => ["face", "car"],
    :pet => ["cat", "dog"],
)


get_sim_formula(sim_type::T) where {T} = get_sim_formula(PredictorStyle(T))

function get_sim_formula(::OnlySplines)
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

function get_sim_formula(::ManyPredictors)
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
function get_sim_formula(::CategoricalPredictors)
    formula = @formula(0 ~ 1 + pet + (1 + pet | subject))
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
    hart = Hartmut()
    source_idx = UnfoldSim.closest_src(hart, "Left Postcentral Gyrus")
    projection = UnfoldSim.magnitude(hart)
    # Only use the first channels/electrodes
    channels = 1:3

    multichannel_signal =
        MultichannelComponent(signal, projection[channels, source_idx], NoNoise())
    return multichannel_signal
end
