import numpy as np
import pandas as pd
import sys

fname = sys.argv[1]
chrom = sys.argv[2]


def GetEnergy(score):
    un_energy = 1.59
    chem_potential = 0
    energy = un_energy+chem_potential + np.log((1-score)/score)
    if energy < 0: return 0
    else: return energy

df = pd.read_csv(fname, sep='\t', header=None, names=['chr','start','end','score'], usecols=[0,1,2,6])
df = df[df["chr"]==chrom]

# Remove duplicates - choose max energy per peak
df = df.groupby(["chr","start","end"], as_index=False).agg{"score": np.max}

# First scale energy to be between 0 and 1
total_reads = np.max(df["score"])
df["score"] = df["score"]/total_reads

# Now convert to energy
df["energy_A"] = df["score"].apply(GetEnergy)

# Set reasonable values for p_ext and p_amp
df['p_ext'] = 0.54 # based on their examples
df['p_amp'] = 0.50 # pcr_rate The mean number of amplified fragments at a location is (1 + p)^n. should be between 0 and 1

df[["chr","start","end","energy_A","p_ext","p_amp"]].to_csv(sys.stdout, index=False, sep="\t")
