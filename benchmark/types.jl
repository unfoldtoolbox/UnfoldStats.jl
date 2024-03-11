abstract type SimulationType end
struct UnitTestSimulation <: SimulationType end
struct BenchmarkSimulation <: SimulationType end

abstract type PredictorStyle end
struct OnlySplines <: PredictorStyle end
struct ManyPredictors <: PredictorStyle end

abstract type ChannelStyle end
struct SingleChannel <: ChannelStyle end
struct MultiChannel <: ChannelStyle end

abstract type EventStyle end
struct SingleEventType <: EventStyle end
struct MultipleEventTypes <: EventStyle end

PredictorStyle(::Type) = OnlySplines()
PredictorStyle(::Type{UnitTestSimulation}) = ManyPredictors()

ChannelStyle(::Type) = SingleChannel()
ChannelStyle(::Type{UnitTestSimulation}) = MultiChannel()

EventStyle(::Type) = SingleEventType()
EventStyle(::Type{UnitTestSimulation}) = MultipleEventTypes()