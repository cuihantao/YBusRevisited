
using YBusRevisited:
    LineUnionParam,
    LineParam,
    LineState,
    LineAlgeb,
    LineRhsF,
    LineRhsG,
    LineJacGy,
    LineConfig
using YBusRevisited: SysConfig, SimConfig

struct LineAux{C}
    expangle::C
    S1_rhs::C
    S2_rhs::C
    yhconj::C
    yhkconj::C
    yhyhkconj::C
    itap2_yhyhkconj::C
    itap_yhkconj::C
end


struct LineRhsGContiguous{F,T}
    rhs::F
    a1_rhs::T
    a2_rhs::T
    v1_rhs::T
    v2_rhs::T
end



struct LineRhsGAddrContiguous{F,T}
    addr::F
    a1_addr::T
    a2_addr::T
    v1_addr::T
    v2_addr::T
end


struct LineJacGyContiguous{F,T}
    mem::F
    _gy1::T
    _gy2::T
    _gy3::T
    _gy4::T
    _gy5::T
    _gy6::T
    _gy7::T
    _gy8::T
    _gy9::T
    _gy10::T
    _gy11::T
    _gy12::T
    _gy13::T
    _gy14::T
    _gy15::T
    _gy16::T
end


struct LineFull{P,S,A,F,GF,GT,GYF,GYT,SV,LA,D,AF,AT}
    n::Int64
    uparam::LineUnionParam
    param::LineParam{P}
    state::LineState{S}
    algeb::LineAlgeb{A}
    rhs_f::LineRhsF{F}
    rhs_g::LineRhsGContiguous{GF,GT}
    jac_gy::LineJacGyContiguous{GYF,GYT}
    service::SV
    aux::LineAux{LA}
    discrete::D
    config::LineConfig
    addr::LineRhsGAddrContiguous{AF,AT}
end

line = sys.data[:Line];

line_gy_store = zeros(16 * line.n);
line_jacgy = LineJacGyContiguous(
    line_gy_store,
    [@view line_gy_store[1+line.n*(i-1):i*line.n] for i = 1:16]...,
);

line_rhs_store = zeros(Float64, 4 * line.n);
line_rhs = LineRhsGContiguous(
    line_rhs_store,
    [@view line_rhs_store[1+line.n*(i-1):i*line.n] for i = 1:4]...,
);

line_addr_store = zeros(Int64, 4 * line.n);
line_addr = LineRhsGAddrContiguous(
    line_addr_store,
    [@view line_addr_store[1+line.n*(i-1):i*line.n] for i = 1:4]...,
);
sys_line_addr = sys.addr[:Line];
line_addr_flatten = collect(
    Iterators.flatten([
        sys_line_addr.algeb.a1_addr,
        sys_line_addr.algeb.a2_addr,
        sys_line_addr.algeb.v1_addr,
        sys_line_addr.algeb.v2_addr,
    ]),
);
line_addr.addr[:] .= line_addr_flatten;


line_aux = LineAux(
    StructArray(zeros(ComplexF64, line.n)),
    StructArray(zeros(ComplexF64, line.n)),
    StructArray(zeros(ComplexF64, line.n)),
    StructArray(zeros(ComplexF64, line.n)),
    StructArray(zeros(ComplexF64, line.n)),
    StructArray(zeros(ComplexF64, line.n)),
    StructArray(zeros(ComplexF64, line.n)),
    StructArray(zeros(ComplexF64, line.n)),
);

lf = LineFull(
    line.n,
    line.uparam,
    line.param,
    line.state,
    line.algeb,
    line.rhs_f,
    line_rhs,
    line_jacgy,
    line.service,
    line_aux,
    line.discrete,
    line.config,
    line_addr,
);

lf.service.bhbhk .= lf.service.bh .+ lf.service.bhk;
lf.service.ghghk .= lf.service.gh .+ lf.service.ghk;

lf.aux.yhconj .= conj.(lf.service.yh);
lf.aux.yhkconj .= conj.(lf.service.yhk);
lf.aux.yhyhkconj .= lf.aux.yhconj .+ lf.aux.yhkconj;

lf.aux.itap2_yhyhkconj .= lf.service.itap2 .* (lf.aux.yhconj .+ lf.aux.yhkconj)
lf.aux.itap_yhkconj .= lf.service.itap .* lf.aux.yhkconj

sys.data[:Line] = lf;
