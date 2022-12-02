# This files contains alternative solutions not used in the paper.

# ------------------------------------------
# Complex number apporoch
# ------------------------------------------
PQYbus = zeros(Float64, 2 * bus_n);
E = exp.(1im * bus_a);
V = E .* bus_v;
S = similar(V);
Ibus = zeros(ComplexF64, bus_n);
buffer = zeros(ComplexF64, bus_n);

@info "Step 1/4: update complex E and V vectors"
b = @benchmark update_E_V($E, $V, $bus_a, $bus_v);
push!(bresults_Ybus, b);

@benchmark begin
    sincos.(bus_a)
end


@info "Step 2/4: calculate line injections"
g_update_Ybus!(PQYbus, S, V, Ybus, bus_n);
b = @benchmark g_update_Ybus!($PQYbus, $S, $V, $Ybus, $bus_n);
push!(bresults_Ybus, b);


# ---------------------------------------------------------------
# A slow version of assembly dS and dR into J
# This method calls `sparse_hcat` and `sparse_vcat` in Julia without observing
# strucural identity

J_dSdR = dSdR_to_J(dS_dVm, dS_dVa);
@benchmark dSdR_to_J($dS_dVm, $dS_dVa)
