using YBusRevisited
using Andes
using BenchmarkTools

Andes.py.config_logger(30)
path = Andes.py.get_case("kundur/kundur_full.xlsx")
ss = Andes.py.run(path, no_output=true, default_config=true)

ss.TDS.init();

sys = YBusRevisited.load_system(ss);

YBusRevisited.g_update!(sys.data[:GENROU], sys.config, sys.sim_config, 0);
YBusRevisited.g_update!(sys.data[:Line], sys.config, sys.sim_config, 0);
