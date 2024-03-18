abstract type SimulationType end

abstract type PredictorStyle end
struct OnlySplines <: PredictorStyle end
struct ManyPredictors <: PredictorStyle end
PredictorStyle(::Type) = OnlySplines()

abstract type ChannelStyle end
struct SingleChannel <: ChannelStyle end
struct MultiChannel <: ChannelStyle end
ChannelStyle(::Type) = SingleChannel()

abstract type EventStyle end
struct SingleEventType <: EventStyle end
struct MultipleEventTypes <: EventStyle end
EventStyle(::Type) = SingleEventType()

struct BenchmarkSimulation <: SimulationType end

struct UnitTestSimulation <: SimulationType end
PredictorStyle(::Type{UnitTestSimulation}) = ManyPredictors()
ChannelStyle(::Type{UnitTestSimulation}) = MultiChannel()
EventStyle(::Type{UnitTestSimulation}) = MultipleEventTypes()

struct AnalysisTutorialSimulation <: SimulationType end
PredictorStyle(::Type{AnalysisTutorialSimulation}) = OnlySplines()
ChannelStyle(::Type{AnalysisTutorialSimulation}) = MultiChannel()
EventStyle(::Type{AnalysisTutorialSimulation}) = SingleEventType()
