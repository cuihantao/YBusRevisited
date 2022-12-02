include("makeYbus.jl");

"""
Complex the complex voltages.
"""
function update_E_V(E, V, bus_a, bus_v)
    @inbounds for m in eachindex(E)
        E[m] = exp(1im * bus_a[m])
        V[m] = E[m] * bus_v[m]
    end
end


function update_E_V(E::StructArray, V, bus_a, bus_v)
    @turbo for m in eachindex(bus_a)
        E.im[m], E.re[m] = sincos(bus_a[m])
        V.re[m], V.im[m] = E.re[m] * bus_v[m], E.im[m] * bus_v[m]
    end
end


function g_update_Ybus!(PQYbus, S, V, Ybus::Union{AbstractSparseMatrix}, bus_n)
    S .= V .* conj(Ybus * V)

    @inbounds for m in eachindex(S)
        PQYbus[m] = real(S[m])
        PQYbus[bus_n+m] = imag(S[m])
    end

    nothing
end


"""
Alternative implementation for network equations using StructArray.

"""
@views function g_update_Ybus_sa!(
    PQYbus,
    Ibus::StructArray,
    S::StructArray,
    Vvec::Vector{ComplexF64},
    V::StructArray,
    Ybus::Union{AbstractSparseMatrix},
    bus_n,
)

    @inbounds for m in eachindex(V)
        Vvec[m] = V[m]
    end

    Ibus .= Ybus * Vvec  # majority of the time; 574 ms on M2

    @turbo Ibus.im .= -Ibus.im
    mul_avx!(S, V, Ibus)

    @turbo PQYbus[1:bus_n] .= S.re
    @turbo PQYbus[bus_n+1:2*bus_n] .= S.im

    nothing
end


"""
g_update using GBMatrix.
"""
function g_update_Ybus_GB!(PQYbus, S, V, Ybus::Union{AbstractSparseMatrix}, bus_n)
    S = V .* conj(Ybus * V)
    PQYbus[1:bus_n] = real(S)
    PQYbus[bus_n+1:2*bus_n] = imag(S)
end
