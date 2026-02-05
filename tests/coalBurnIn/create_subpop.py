#%%
import msprime as msp
import pyslim as pys
from IPython.display import SVG, display 

demography = msp.Demography()
demography.add_population(name="p12", initial_size=100)
demography.add_population(name="p34", initial_size=100)
demography.add_population(name="p56", initial_size=100)
demography.add_population(name="p3456", initial_size=100)
demography.add_population(name="pbase", initial_size=100)

demography.add_population(name="p1", initial_size=100)
demography.add_population(name="p2", initial_size=100)
demography.add_population(name="p3", initial_size=100)
demography.add_population(name="p4", initial_size=100)
demography.add_population(name="p5", initial_size=100)
demography.add_population(name="p6", initial_size=100)

demography.add_population_split(time=1000, derived = ["p5", "p6"], ancestral="p56")
demography.add_population_split(time=1500, derived = ["p3", "p4"], ancestral="p34")
demography.add_population_split(time=2000, derived = ["p34", "p56"], ancestral="p3456")
demography.add_population_split(time=2500, derived = ["p1", "p2"], ancestral="p12")
demography.add_population_split(time=3000, derived = ["p12", "p3456"], ancestral="pbase")

breaks = [0, 3_000_000, 6_000_000, 9_000_000]
recomb_map = msp.RateMap(
    position = breaks,
    rate = [5e-8, 0.5e-8, 5e-8]
)

ts = msp.sim_ancestry(samples={"p1": 5000, "p2": 5000, "p3": 5000,
                               "p4": 5000, "p5": 5000, "p6": 5000}, 
                      demography = demography,
                      recombination_rate=recomb_map)
#ts
#SVG(ts.draw_svg(y_axis=True))
ots = pys.annotate(ts, model_type="WF", tick=1, stage="late")

mut_map = msp.RateMap(
    position = breaks,
    rate = [1e-8, 1e-9, 1e-8]
)

ots = msp.sim_mutations(ots, rate = mut_map, 
                        model = msp.SLiMMutationModel(type = 1),
                        keep = True
                  )
#SVG(ots.draw_svg(y_axis=True))
pts = ots.tables.tree_sequence()
pts.dump("/mnt/c/GitHub/SLiMTests/tests/coalBurnIn/test.trees")
# %%
