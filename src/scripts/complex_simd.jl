using MKL
using LoopVectorization, LinearAlgebra, StructArrays, BenchmarkTools, Test
using LoopVectorization: vpermilps177, vmovshdup, vfmsubadd, vfmaddsub, vmovsldup


# BLAS.set_num_threads(1); @show BLAS.get_config()

const MatrixFInt64 = Union{SubArray,Vector{Float64},Vector{Int}}

@inline function mul_avx!(C::MatrixFInt64, A::MatrixFInt64, B::MatrixFInt64)
    @turbo for m in eachindex(A)
        C[m] = A[m] * B[m]
    end
end

@inline function mul_add_avx!(C::MatrixFInt64, A::MatrixFInt64, B::MatrixFInt64, factor = 1)
    @turbo for m in eachindex(A)
        C[m] += factor * A[m] * B[m]
    end
end


const StructMatrixComplexFInt64 =
    Union{StructArray{ComplexF64,2},StructArray{Complex{Int},2}}

function mul_avx!(C, A, B)
    mul_avx!(C.re, A.re, B.re)     # C.re = A.re * B.re
    mul_add_avx!(C.re, A.im, B.im, -1) # C.re = C.re - A.im * B.im
    mul_avx!(C.im, A.re, B.im)     # C.im = A.re * B.im
    mul_add_avx!(C.im, A.im, B.re)     # C.im = C.im + A.im * B.re
end
