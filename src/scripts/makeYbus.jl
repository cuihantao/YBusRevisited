"""
This script contains functions to make Ybus admittance matrix.

At the end, this script creates an admittance matrix in the variable `Ybus`.
"""

using YBusRevisited: LineRhsGAddr
using SparseArrays
using LinearAlgebra

include("linefull.jl");


"""
Build admittance matrix for System (Line and Bus).
"""
function Ymatrix(line::LineFull, line_addr::LineRhsGAddr{T}, n_bus) where {T}
    rows = zeros(Int64, 4line.n)
    cols = zeros(Int64, 4line.n)
    vals = zeros(Complex{Float64}, 4line.n)

    (; u, tap, phi) = line.param
    (; gh, bh, gk, bk, ghk, bhk, itap, itap2) = line.service

    (; a1_addr, a2_addr) = line_addr

    y1 = u .* (gh + bh * 1im)
    y2 = u .* (gk + bk * 1im)
    y12 = u .* (ghk + bhk * 1im)

    tapc = tap .* exp.(1im * phi)
    itapc = 1 ./ tapc
    itapcconj = conj.(itapc)

    Ymatrix!(
        line.n,
        a1_addr,
        a2_addr,
        rows,
        cols,
        vals,
        y1,
        y2,
        y12,
        itapc,
        itap2,
        itapcconj,
    )

    sparse(rows, cols, vals, n_bus, n_bus)
end


"""
Allocation-free function for Ymatrix building
"""
Base.@inline function Ymatrix!(
    n_line,
    a1_addr,
    a2_addr,
    rows::Vector{Int64},
    cols::Vector{Int64},
    vals::Vector{Complex{Float64}},
    y1,
    y2,
    y12,
    itapc,  # itap * e^(-1im * phi)
    itap2,
    itapcconj,
)

    @inbounds for i in eachindex(a1_addr)

        rows[i] = a1_addr[i]
        cols[i] = a1_addr[i]
        vals[i] = (y1[i] + y12[i]) * itap2[i]

        rows[n_line+i] = a1_addr[i]
        cols[n_line+i] = a2_addr[i]
        vals[n_line+i] = -y12[i] * itapcconj[i]

        rows[2n_line+i] = a2_addr[i]
        cols[2n_line+i] = a1_addr[i]
        vals[2n_line+i] = -y12[i] * itapc[i]

        rows[3n_line+i] = a2_addr[i]
        cols[3n_line+i] = a2_addr[i]
        vals[3n_line+i] = y12[i] + y2[i]

    end
    nothing
end

bus = sys.data[:Bus];

Ybus = Ymatrix(lf, sys.addr[:Line].rhs_g, bus.n);

bus_n = bus.n;
bus_v = get_tmp(bus.algeb.v, 0);
bus_a = get_tmp(bus.algeb.a, 0);

U = exp.(1im * bus_a);
V = U .* bus_v;
Ic = Ybus * V;

diagVc = spdiagm(bus_n, bus_n, V);
diagVn = spdiagm(bus_n, bus_n, U);
diagIc = spdiagm(bus_n, bus_n, Ic);
