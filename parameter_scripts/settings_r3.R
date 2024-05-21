# Setting initial values
ini <- list(n_bac = 100,                 # initial number of bacteria
            max_bac = 10000,             # max bacteria in the simulation
            n_grid = 100,                # number of cells in x, y in spatial grid
            max_bac_cell = 50,            # max number of bacteria allowed per cell
            resistance_perc = 0,         # percent resistant bacteria at start
            antibiotic_conc = 0)         # initial concentration of antibiotics

# Setting parameters
p <- list(conjugation_prob = 0.05,       # probability of successful HGT
          protective_factor = 1,         # number of bacteria protected by each resistant bacterium
          mut_ben_prob = 1e-4,           # 10^-8 mutation rate for beneficial mutations for non-mutators
          mut_del_prob = 1,              # 10^-4 mutation rate for deleterious mutations for non-mutators
          mut_lethal_prob = 0.1,         # 10^-5 mutation rate for lethal mutations for non-mutators
          mut_rate = 0.005,              # based on 10^-9 mutations per base in E. coli and a 5 Mb genome
          
          # Antibiotic parameters
          antibiotic_MIC = 4,            # MIC of the antibiotic for the bacteria in question
          dose_conc = 0,                 # increase in concentration of antibiotics after a dose
          dose_first = con_factor*0,     # time at which the first dose is given
          dose_interval = 0,             # time interval (hours) between doses
          dose_number = 1,               # number of doses given
          half_life = con_factor*1)      # half life of the antibiotic in the system
