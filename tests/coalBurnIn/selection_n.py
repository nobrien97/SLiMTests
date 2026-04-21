#%%
#from Nick
import msprime as msp
import pyslim as pys
from IPython.display import SVG, display 

dmg = msp.Demography()
dmg.add_population(name="ancestral", initial_size=50000)
dmg.add_population(name="pop1", initial_size=10000)
dmg.add_population(name="pop2", initial_size=10000)
dmg.add_population(name="pop3", initial_size=10000)
dmg.add_population(name="pop4", initial_size=10000)
dmg.add_population(name="pop5", initial_size=10000)
dmg.add_population(name="pop6", initial_size=10000)
dmg.add_population(name="pop7", initial_size=10000)
dmg.add_population(name="pop8", initial_size=10000)

dmg.add_population_split(time=3000, derived=["pop3", "pop4"], ancestral="pop1")
dmg.add_population_split(time=2000, derived=["pop5", "pop6"], ancestral="pop2")
dmg.add_population_split(time=1000, derived=["pop7", "pop8"], ancestral="pop6")

dmg.add_population_split(
    time=5000,
    derived=["pop1", "pop2"],
    ancestral="ancestral"
)

dmg.sort_events()

#recombination rates
LOW_REC_REGIONS = [(250_000, 300_000), (370_000, 380_000),(600_000, 630_000), (700_000, 800_000)]
positions = [0]
rates = []

for start, end in LOW_REC_REGIONS:
    positions.extend([start, end])
    rates.extend([1e-8, 1e-9])

positions.append(2000000)
rates.append(1e-8)

recomb_map = msp.RateMap(position=positions, rate=rates)



# breaks = [0, 1_000_000, 1_500_000, 2_000_000]
# recomb_map = msp.RateMap(
#     position = breaks,
#     rate = [5e-8, 0.5e-8, 5e-8]
# )
N_SAMPLES=80
ts = msp.sim_ancestry(
    samples={"pop3": N_SAMPLES, "pop4": N_SAMPLES, "pop5": N_SAMPLES, "pop7": N_SAMPLES, "pop8": N_SAMPLES},
    demography=dmg,
    recombination_rate=recomb_map,
    sequence_length=2000000,)


ts = pys.annotate(ts, model_type="WF", tick=1, stage="late")

# mut_map = msp.RateMap(
#     position = positions ,
#     rate = 1e-8
# )


ts = msp.sim_mutations(ts, rate = 1e-8 , 
                        model = msp.SLiMMutationModel(type = 1),
                        keep = True
                  )
ts.dump("new_chapt4.trees")
#%%