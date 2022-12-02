# ------------------
# Test cases
# ------------------


if @isdefined(case) == false || case == ""
    case = "$(HOME_DIR)/repos/matpower/data/case9241pegase.m"
    @info "`case` not defined. Using default case9241pegase"
end

ss = Andes.py.load(case, no_output = true, default_config = true);

ss.PFlow.init();

sys = YBusRevisited.load_system(ss);
xy0 = ss.dae.xy;
out = similar(ss.dae.xy);

res!(out, xy0)
