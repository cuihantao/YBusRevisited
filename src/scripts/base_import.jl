using Revise
using YBusRevisited

using BenchmarkTools
using LoopVectorization
using Andes
using PreallocationTools
using StructArrays

using LinearAlgebra
using MKL

BLAS.set_num_threads(1);

HOME_DIR = ENV["HOME"];

# ------------------
#  Residual methods
# ------------------

res! = (output, xy) -> YBusRevisited.evaluate_residuals!(sys, output, xy);
