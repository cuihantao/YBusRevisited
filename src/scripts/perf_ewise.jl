using LinuxPerf: make_bench, enable!, disable!, reset!, reasonable_defaults, counters

const bench = make_bench(reasonable_defaults);

@noinline function perf_gy_update_optimized2!(bench, lf, sys_config, sys_sim_config)
    enable!(bench)
    for i = 1:1000

        gy_update_turbo!(lf, sys_config, sys_sim_config, 0)

    end
    disable!(bench)
end

reset!(bench)
perf_gy_update_optimized2!(bench, lf, sys.config, sys.sim_config)

counters(bench)

#=

Linux Perf of gy on 12900k ; case_SyntheticUSA

┌───────────────────────┬───────────────┬─────────────┐
│                       │ Events        │ Active Time │
├───────────────────────┼───────────────┼─────────────┤
│             hw:cycles │ 4,577,698,819 │ 100.0 %     │
├───────────────────────┼───────────────┼─────────────┤
│       hw:cache_access │ 581,952,974   │ 100.0 %     │
├───────────────────────┼───────────────┼─────────────┤
│       hw:cache_misses │ 60,248,850    │ 100.0 %     │
├───────────────────────┼───────────────┼─────────────┤
│           hw:branches │ 19,537,258    │ 100.0 %     │
├───────────────────────┼───────────────┼─────────────┤
│ hw:branch_mispredicts │ 2,538         │ 100.0 %     │
├───────────────────────┼───────────────┼─────────────┤
│       hw:instructions │ 3,690,449,376 │ 100.0 %     │
├───────────────────────┼───────────────┼─────────────┤
│       sw:ctx_switches │ 0             │ 100.0 %     │
├───────────────────────┼───────────────┼─────────────┤
│        sw:page_faults │ 0             │ 100.0 %     │
├───────────────────────┼───────────────┼─────────────┤
│  sw:minor_page_faults │ 0             │ 100.0 %     │
├───────────────────────┼───────────────┼─────────────┤
│  sw:major_page_faults │ 0             │ 100.0 %     │
├───────────────────────┼───────────────┼─────────────┤
│     sw:cpu_migrations │ 0             │ 100.0 %     │
└───────────────────────┴───────────────┴─────────────┘

=#
