include("g_Ybus.jl");
include("gy_csc.jl");
include("gy_elemwise.jl");
include("makeYbus.jl");
include("complex_simd.jl")


@info "Benchmark: Ybus method for $case"

bresults_Ybus = Vector{BenchmarkTools.Trial}();

# -----------------------------------------------------
# Residuals (power injection) calculation
#


# --- Alternative StructArray storage for complex numbers ---
# NOTE: this approach is inefficient for Y * V because of memory layout

E = exp.(1im * bus_a);

PQYbus = (zeros(Float64, 2 * bus_n));
V = StructArray(V);
E = StructArray(E);
S = similar(V);
Ibus = StructArray(zeros(ComplexF64, bus_n));
buffer = StructArray(zeros(ComplexF64, bus_n));

@info "Step 1/4: update complex E and V vectors"
update_E_V(E, V, bus_a, bus_v);
b = @benchmark update_E_V($E, $V, $bus_a, $bus_v);
push!(bresults_Ybus, b);

# ---- Alternative version to calculate network inj ---

# NOTE: must use `Ybus * Vvec` becase `Ybus * V` is extremely slow
Vvec = zeros(ComplexF64, bus_n);
@info "Step 2/4: calculate line injections"
g_update_Ybus_sa!(PQYbus, S, Ibus, Vvec, V, Ybus, bus_n);
b = @benchmark g_update_Ybus_sa!($PQYbus, $S, $Ibus, $Vvec, $V, $Ybus, $bus_n);
push!(bresults_Ybus, b);

# ----------------------------------------------
# ---- Jacobian part using CSC matrix----

Yx = Ybus.nzval;
Yp = Ybus.colptr;
Yi = Ybus.rowval;
bus_n = length(V);

dS_dVm = copy(Ybus);
dS_dVa = copy(Ybus);

dS_dVm, dS_dVa = dSbus_dV_csc(Ybus, dS_dVm, dS_dVa, buffer, Ibus, Yx, Yp, Yi, V, E);

@info "Step 3/4: calculate dS_dVm and dS_dVa"
b = @benchmark dS_dVm, dS_dVa =
    dSbus_dV_csc($Ybus, $dS_dVm, $dS_dVa, $buffer, $Ibus, $Yx, $Yp, $Yi, $V, $E);
push!(bresults_Ybus, b);

# -------------------------------------------------
#
# Building sparse Jacobian matrix
#
# -------------------------------------------------
# Fast assembly of sparse matrix

nnz_Ybus = nnz(Ybus);
J_dSdR = dSdR_to_J_keep_zeros(dS_dVm, dS_dVa); # drops zero elements
Jshape = copy(J_dSdR);

# J1
perm_J1 = zeros(Int64, nnz_Ybus);
rowidx, colidx, nzvals = findnz(dS_dVa);
find_nzloc!(perm_J1, Jshape, rowidx, colidx);

# J2
perm_J2 = copy(perm_J1);
inc_perm_idx_col!(perm_J2, colidx, Ybus);

# J3
perm_J3 = copy(perm_J1);
perm_J3 .+= 2 * nnz_Ybus;

# J4
perm_J4 = copy(perm_J3);
inc_perm_idx_col!(perm_J4, colidx, Ybus);

@info "Step 4/4: Assemble dS_dVm and dS_dVa into J"

assemble_J_Ybus(Jshape, dS_dVa, dS_dVm, perm_J1, perm_J2, perm_J3, perm_J4)
b = @benchmark assemble_J_Ybus(Jshape, dS_dVa, dS_dVm, perm_J1, perm_J2, perm_J3, perm_J4);
push!(bresults_Ybus, b);

@show bresults_Ybus
