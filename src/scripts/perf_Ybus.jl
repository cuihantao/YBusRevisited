# -------------------------------------------------------
# LinuxPerf benchmark
# -------------------------------------------------------

using LinuxPerf: make_bench, enable!, disable!, reset!, reasonable_defaults, counters

const bench = make_bench(reasonable_defaults);

# b = @benchmark dS_dVm, dS_dVa = dSbus_dV_csc($Ybus, $dS_dVm, $dS_dVa, $buffer, $Ibus, $Yx, $Yp, $Yi, $V, $E);

@noinline function perf_dSbus_dV_csc(
    bench,
    Ybus,
    dS_dVm,
    dS_dVa,
    buffer,
    Ibus,
    Yx,
    Yp,
    Yi,
    V,
    E,
)
    enable!(bench)
    for i = 1:1000
        dSbus_dV_csc(Ybus, dS_dVm, dS_dVa, buffer, Ibus, Yx, Yp, Yi, V, E)
    end
    disable!(bench)
end

reset!(bench)
perf_dSbus_dV_csc(bench, Ybus, dS_dVm, dS_dVa, buffer, Ibus, Yx, Yp, Yi, V, E);

counters(bench)

#=

Perf data on 12900k: case_SyntheticUSA with 1000 runs

┌───────────────────────┬────────────────┬─────────────┐
│                       │ Events         │ Active Time │
├───────────────────────┼────────────────┼─────────────┤
│             hw:cycles │ 14,327,078,176 │ 100.0 %     │
├───────────────────────┼────────────────┼─────────────┤
│       hw:cache_access │ 828,580,410    │ 100.0 %     │
├───────────────────────┼────────────────┼─────────────┤
│       hw:cache_misses │ 156,647,072    │ 100.0 %     │
├───────────────────────┼────────────────┼─────────────┤
│           hw:branches │ 5,474,516,729  │ 100.0 %     │
├───────────────────────┼────────────────┼─────────────┤
│ hw:branch_mispredicts │ 173,767,459    │ 100.0 %     │
├───────────────────────┼────────────────┼─────────────┤
│       hw:instructions │ 50,308,631,851 │ 100.0 %     │
├───────────────────────┼────────────────┼─────────────┤
│       sw:ctx_switches │ 0              │ 100.0 %     │
├───────────────────────┼────────────────┼─────────────┤
│        sw:page_faults │ 0              │ 100.0 %     │
├───────────────────────┼────────────────┼─────────────┤
│  sw:minor_page_faults │ 0              │ 100.0 %     │
├───────────────────────┼────────────────┼─────────────┤
│  sw:major_page_faults │ 0              │ 100.0 %     │
├───────────────────────┼────────────────┼─────────────┤
│     sw:cpu_migrations │ 0              │ 100.0 %     │
└───────────────────────┴────────────────┴─────────────┘

=#
