import numpy as np
import pandas as pd
import sys

fname = sys.argv[1]
pcr_rate = float(sys.argv[2])

un_energy = 1.59
chem_potential = 0

df = pd.read_csv(fname, sep='\t', header=None, names=['chr','start','end','energy_A'], usecols=[0,1,2,6])

#total_reads = np.sum(df['energy_A'].tolist())
total_reads = df['energy_A'][0] +1
df['energy_A'] = df['energy_A']/total_reads
df['energy_A'] = [ i if i > 0 else 0 for i in un_energy + chem_potential + np.log(1-df['energy_A'])/df['energy_A']]
df['p_ext'] = 0.54
df['p_amp'] = pcr_rate

df.to_csv(sys.stdout, index=False, sep="\t")
