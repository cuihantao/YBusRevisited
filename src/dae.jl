"""
Numerical DAE structs.
"""
struct DAE{T}
    m::Int64
    n::Int64

    x::T  # dualcache -> DiffCache     get_tmp(x, NewValue)
    y::T

    f::T
    g::T
end
