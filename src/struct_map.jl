include("models/models.jl")

model_map = Dict{Symbol,Type}(
    :Bus => Bus,
    :PQ => PQ,
    :PV => PV,
    :Slack => Slack,
    :Shunt => Shunt,
    :Line => Line,
)

fields_map = Dict{Symbol,Vector}(
    :Bus => [
        BusUnionParam,
        BusParam,
        BusState,
        BusAlgeb,
        BusRhsF,
        BusRhsG,
        BusService,
        BusDiscrete,
        BusConfig,
    ],
    :PQ => [
        PQUnionParam,
        PQParam,
        PQState,
        PQAlgeb,
        PQRhsF,
        PQRhsG,
        PQService,
        PQDiscrete,
        PQConfig,
    ],
    :PV => [
        PVUnionParam,
        PVParam,
        PVState,
        PVAlgeb,
        PVRhsF,
        PVRhsG,
        PVService,
        PVDiscrete,
        PVConfig,
    ],
    :Slack => [
        SlackUnionParam,
        SlackParam,
        SlackState,
        SlackAlgeb,
        SlackRhsF,
        SlackRhsG,
        SlackService,
        SlackDiscrete,
        SlackConfig,
    ],
    :Shunt => [
        ShuntUnionParam,
        ShuntParam,
        ShuntState,
        ShuntAlgeb,
        ShuntRhsF,
        ShuntRhsG,
        ShuntService,
        ShuntDiscrete,
        ShuntConfig,
    ],
    :Line => [
        LineUnionParam,
        LineParam,
        LineState,
        LineAlgeb,
        LineRhsF,
        LineRhsG,
        LineService,
        LineDiscrete,
        LineConfig,
    ],
)

address_map = Dict{Symbol,Type}(
    :Bus => BusAddr,
    :PQ => PQAddr,
    :PV => PVAddr,
    :Slack => SlackAddr,
    :Shunt => ShuntAddr,
    :Line => LineAddr,
)

address_substructs = Dict{Symbol,Vector}(
    :Bus => [BusStateAddr, BusAlgebAddr, BusRhsFAddr, BusRhsGAddr],
    :PQ => [PQStateAddr, PQAlgebAddr, PQRhsFAddr, PQRhsGAddr],
    :PV => [PVStateAddr, PVAlgebAddr, PVRhsFAddr, PVRhsGAddr],
    :Slack => [SlackStateAddr, SlackAlgebAddr, SlackRhsFAddr, SlackRhsGAddr],
    :Shunt => [ShuntStateAddr, ShuntAlgebAddr, ShuntRhsFAddr, ShuntRhsGAddr],
    :Line => [LineStateAddr, LineAlgebAddr, LineRhsFAddr, LineRhsGAddr],
)
