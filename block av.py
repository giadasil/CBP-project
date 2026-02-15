import MDAnalysis as mda
from MDAnalysis.analysis import distances
import numpy as np

def get_block_sem(data, block_size=1000):
    """Calculates mean and SEM using block averaging."""
    if len(data) == 0: return 0, 0
    n_blocks = len(data) // block_size
    if n_blocks == 0: return np.mean(data), 0
    
    # Reshape into blocks and calculate the mean of each block
    block_means = data[:n_blocks*block_size].reshape(n_blocks, block_size).mean(axis=1)
    return np.mean(data), np.std(block_means) / np.sqrt(n_blocks)

def calculate_distances(gro, xtc):
    """Extracts distances for Lys27-His102 CA atoms over time."""
    u = mda.Universe(gro, xtc)
    sel1 = u.select_atoms("resid 27 and name CA")
    sel2 = u.select_atoms("resid 102 and name CA")
    
    dist_list = []
    time_list = []
    
    for ts in u.trajectory:
        # Calculate distance between the two atom groups
        d = distances.dist(sel1, sel2)[2][0]
        dist_list.append(d)
        time_list.append(ts.time / 1000.0) # Convert ps to ns
        
    return np.array(dist_list), np.array(time_list)

# --- Process WT and AF (Full Trajectory) ---
systems = {
    "Wild Type": ("wt_protein_only.gro", "wt_protein_only.xtc"),
    "AF Mutant": ("af_protein_only.gro", "Alpha_protein_only.xtc")
}

print(f"{'System':<20} | {'Mean (Å)':<10} | {'SEM':<8}")
print("-" * 45)

for name, files in systems.items():
    dists, _ = calculate_distances(files[0], files[1])
    avg, sem = get_block_sem(dists)
    print(f"{name:<20} | {avg:<10.3f} | {sem:<8.3f}")

# --- Process Chimera (Split at 350 ns) ---
c_dist, c_time = calculate_distances("chimera_protein_only.gro", "chimera_protein_only.xtc")

# Boolean masks for pre and post 350ns
pre_mask = c_time <= 350
post_mask = c_time > 350

avg_pre, sem_pre = get_block_sem(c_dist[pre_mask])
avg_post, sem_post = get_block_sem(c_dist[post_mask])

print(f"{'Chimera (Pre-350ns)':<20} | {avg_pre:<10.3f} | {sem_pre:<8.3f}")
print(f"{'Chimera (Post-350ns)':<20} | {avg_post:<10.3f} | {sem_post:<8.3f}")

import matplotlib.pyplot as plt
import seaborn as sns

def plot_results(times, dists, name):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Time Series
    ax1.plot(times, dists, alpha=0.5, color='royalblue')
    ax1.axvline(x=350, color='red', linestyle='--', label='Jump Point (350ns)')
    ax1.set_xlabel('Time (ns)')
    ax1.set_ylabel('Distance (Å)')
    ax1.set_title(f'{name} Distance over Time')
    ax1.legend()

    # Distribution (Pre vs Post)
    pre = dists[times <= 350]
    post = dists[times > 350]
    sns.kdeplot(pre, ax=ax2, fill=True, label='Pre-Jump', color='green')
    sns.kdeplot(post, ax=ax2, fill=True, label='Post-Jump', color='orange')
    ax2.set_xlabel('Distance (Å)')
    ax2.set_title('Distance Distribution Population')
    ax2.legend()

    plt.tight_layout()
    plt.show()

# Run this using the c_dist and c_time variables from the previous script
plot_results(c_time, c_dist, "Chimera")