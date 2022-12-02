# -------------------------------------------------
# Test power injections match
#
# -------------------------------------------------

# Test power injections match
@assert (maximum(abs.(PQYbus - PQelem)) < 1e-10) "Calculated power injections do not match"

# -------------------------------------------------
# Test Jacobians match
# -------------------------------------------------

# --- ground truth from matrix multiplication---
# Same or similar equations can be found in PSAT

dS = diagVc * conj(Ybus * diagVn) .+ conj(diagIc) * diagVn;    # dS_dVm
dR = 1im * conj(conj(diagVc) * (diagIc - Ybus * diagVc));      # dS_dVa

# @benchmark begin
#     $dS = $diagVc * conj($Ybus * $diagVn) .+ conj($diagIc) * $diagVn;    # dS_dVm
#     $dR = 1im * conj(conj($diagVc) * ($diagIc - $Ybus * $diagVc));      # dS_dVa
# end

# V = Vector{ComplexF64}(V);
# Ic = Vector{ComplexF64}(Ic);
# dR = 1im * diagVc * (conj(diagIc) - conj(Ybus * diagVc));

@assert nnz(droptol!(dS_dVm - dS, 1e-8)) == 0 "dS_dVm and dS do not match"
@assert nnz(droptol!(dS_dVa - dR, 1e-8)) == 0 "dS_dVa and dR do not match"

@assert nnz(droptol!(Jshape - Jelem, 1e-8)) == 0 "Jshape from dS_dVa and dS_dVm does not match Jelem"
@assert (nnz(droptol!(Jelem - J_dSdR, 1e-8)) == 0) "Jelem does not match J_dSdR from matrix calc."
@assert (nnz(droptol!(Jshape - J_dSdR, 1e-8)) == 0) "Jshape does not match J_dSdR from matrix calc."
