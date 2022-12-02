
# --- separate a1_addr and `a2_addr` into unique and nonunique indices into
# PQinj ---
# probably not necessary, because it's fast enough

function find_perm_g!(l_nonunique, r_nonunique, rperm, addr)
    rperm .= 0
    for m = 1:length(rperm)
        rperm[m] = findfirst(isequal(m), addr)
    end
    r_nonunique = setdiff(1:length(addr), rperm)
    l_nonunique = addr[r_nonunique]

    nothing
end


sys_line_addr = sys.addr[:Line];
PQ_eqn_addr = collect(
    Iterators.flatten([
        sys_line_addr.algeb.a1_addr,
        sys_line_addr.algeb.a2_addr,
        sys_line_addr.algeb.v1_addr,
        sys_line_addr.algeb.v2_addr,
    ]),
);


perm_line_rhs = zeros(Int64, 2 * bus_n);
PQ_eqn_nonunique = zeros(Int64, 4 * lf.n - 2 * bus_n);  # goes to left-hand side, i.e., PQinj[eqn_nonunique] = ...
line_rhs_nonunique = similar(PQ_eqn_nonunique);          # use to index into [a1_rhs, a2_rhs, ...]
find_perm_g!(PQ_eqn_nonunique, line_rhs_nonunique, perm_line_rhs, PQ_eqn_addr);



# -------------------- Another version is here. Tested slower ----

# g_perm = zeros(Int64, 2*bus_n);
# gl_nonunique = zeros(Int64, 4*line.n - 2*bus_n);
# gr_nonunique = zeros(Int64, 4*line.n - 2*bus_n);

# find_perm_g!(gl_nonunique, gr_nonunique, g_perm, lf.addr.addr)
# call_collect_inj(PQelem, lf, g_perm, gl_nonunique, gr_nonunique);
# b = @benchmark call_collect_inj($PQelem, $lf, $g_perm, $gl_nonunique, $gr_nonunique);


# --- Assemble `gy` by indexing into `Jelem.nzval` ---
# Turns out slower than the two-step method.

perm = zeros(Int64, length(gy_rows));
find_nzloc!(perm, Jelem, gy_rows, gy_cols)


assemble_J_elem!(Jelem, gy_vals, perm)
# @benchmark assemble_J_elem!($Jelem, $gy_vals, $perm)


Jelem = Jelem[1:25_000, 1:25_000];
PQelem = PQelem[1:25_000];


using Random
using KLU
factor = klu(Jelem);
klu!(factor, Jelem);
Jelem \ PQelem;

@btime klu($Jelem);
@btime klu!($factor, $Jelem);

q = repeat(PQelem, 1, 20000);

retained = rand((1:20_000), 50_000);
qz = zeros(size(q));

for (idx, col) in enumerate(retained)
    row = Int(ceil(idx / 2))
    qz[row, col] = q[row, col]
end

qzs = sparse(qz);
@benchmark begin
    $factor \ $qzs
end


using Pardiso

ps = MKLPardisoSolver()

@benchmark solve($ps, $Jelem, $PQelem)

factor = klu(Jelem);
factor \ PQelem

Jelem \ PQelem

@benchmark begin
    klu($Jelem) \ $PQelem
end
