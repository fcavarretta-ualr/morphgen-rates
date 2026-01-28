import pandas as pd
from morphgen_rates import compute_rates

# Load Sholl plot summary statistics (bin counts + variance) from CSV
df_sholl = pd.read_csv("data/pyr_apical_sholl_plot.csv", index_col=0)

# manipulate the data
df_sholl = df_sholl.T.describe().T[['mean', 'std']]
bin_size = df_sholl.index[1] - df_sholl.index[0]
df_sholl = df_sholl.to_numpy()

# Load bifurcation summary statistics from CSV
df_bif = pd.read_csv("data/pyr_bifurcations.csv", index_col=0).to_numpy()


# Bundle inputs exactly as loaded (no preprocessing)
data = {
  "sholl": {
    'bin_size':bin_size, 
    'mean':df_sholl[:-1, 0],
    'var':df_sholl[:-1, 1] ** 2,
    },
  "bifurcations": {
    'mean':df_bif.mean(),
    'var':df_bif.var()
    },
}

max_step_size = 5.

rates = compute_rates(data, max_step_size=max_step_size)

print("Bifurcation rate:", rates.get("bifurcation_rate"))
print("Annihilation rate:", rates.get("annihilation_rate"))
