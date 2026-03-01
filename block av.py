import MDAnalysis as mda
from MDAnalysis.analysis import distances
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

u_wt = mda.Universe("wt_protein_only.gro", "wt_protein_only.xtc")
u_af = mda.Universe("af_protein_only.gro", "Alpha_protein_only.xtc")
u_ch = mda.Universe("chimera_protein_only.gro", "chimera_protein_only.xtc")

total_frames = len(u_wt.trajectory)
start_frame = total_frames - 9000

def get_distance_array(u, start, res1, res2):
    # Extracts Ca-Ca distance for every frame in the slice
    sel1 = u.select_atoms(f"resid {res1} and name CA")
    sel2 = u.select_atoms(f"resid {res2} and name CA")
    
    distances = []
    for ts in u.trajectory[start:]:
        # Calculate Euclidean distance between the two atom positions
        dist = np.linalg.norm(sel1.positions - sel2.positions)
        distances.append(dist)
    return np.array(distances)

# Extract the full 9000-frame arrays
# start_frame was defined as total_frames - 9000
dist_wt_9000 = get_distance_array(u_wt, start_frame, 102, 104)
dist_af_9000 = get_distance_array(u_af, start_frame, 102, 104)
dist_ch_9000 = get_distance_array(u_ch, start_frame, 102, 104)

# Pass these arrays into block stats function
def get_block_stats(data, num_blocks=10):
    block_size = len(data) // num_blocks
    trimmed_data = data[:block_size * num_blocks]
    blocks = trimmed_data.reshape(num_blocks, block_size)
    block_means = np.mean(blocks, axis=1) # mean for wach block
    return np.mean(block_means), np.std(block_means) / np.sqrt(num_blocks) # total mean and SEM

# Calculate final stats for the bar chart
mean_wt, err_wt = get_block_stats(dist_wt_9000)
mean_af, err_af = get_block_stats(dist_af_9000)
mean_ch, err_ch = get_block_stats(dist_ch_9000)

print(f"Wild Type: {mean_wt:.3f} ± {err_wt:.3f} Å")
print(f"Chimera:   {mean_ch:.3f} ± {err_ch:.3f} Å")
print(f"AF Mutant: {mean_af:.3f} ± {err_af:.3f} Å")

systems = ['Wild Type', 'AlphaFold Mutant', 'Chimera Mutant']
data_list = [dist_wt_9000, dist_af_9000, dist_ch_9000]
colors = ['#1f77b4', '#2ca02c', '#d62728'] 
means = []
errors = []

for d in data_list:
    m, e = get_block_stats(d, num_blocks=10)
    means.append(m)
    errors.append(e)

plt.figure(figsize=(8, 6))
bars = plt.bar(systems, means, yerr=errors, capsize=10, color=colors, alpha=0.7, edgecolor='black')

plt.ylabel('Distance ($C\\alpha$-$C\\alpha$) ($\AA$)', fontsize=12)
plt.title('Block-Averaged Active Site Distances (Res 102-104)', fontsize=14)
plt.ylim(min(means) - 0.5, max(means) + 0.5) # Zoom in on the differences
plt.grid(axis='y', linestyle='--', alpha=0.3)

for bar in bars:
    yval = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2, yval + 0.05, f'{yval:.2f} $\AA$', ha='center', va='bottom')

plt.tight_layout()
plt.show()

# Generate the 10 independent means for each system
def get_block_means(data, num_blocks=10):
    block_size = len(data) // num_blocks
    blocks = data[:block_size * num_blocks].reshape(num_blocks, block_size)
    return np.mean(blocks, axis=1)
blocks_wt = get_block_means(dist_wt_9000)
blocks_af = get_block_means(dist_af_9000)
blocks_ch = get_block_means(dist_ch_9000)
# Welch's t-test for statistical significance (equal_var=False)
t_stat_ch, p_val_ch = stats.ttest_ind(blocks_wt, blocks_ch, equal_var=False)
t_stat_af, p_val_af = stats.ttest_ind(blocks_wt, blocks_af, equal_var=False)

print(f"--- Welch's T-test Results (Comparison to Wild Type) ---")
print(f"Chimera Mutant:  t = {t_stat_ch:.3f},  p = {p_val_ch:.2e}")
print(f"AF Mutant:       t = {t_stat_af:.3f},  p = {p_val_af:.2e}")

# Scientific interpretation
def interpret_p(p):
    if p < 0.001: return "*** Highly Significant"
    if p < 0.01:  return "** Very Significant"
    if p < 0.05:  return "* Significant"
    return "Not Significant"

print(f"\nChimera Significance: {interpret_p(p_val_ch)}")
print(f"AF Mutant Significance: {interpret_p(p_val_af)}")
