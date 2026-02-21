import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

def calculate_rg(gro, xtc):
    u = mda.Universe(gro, xtc)
    protein = u.select_atoms("name CA")
    rg_values = []
    
    for ts in u.trajectory:
        rg_values.append(protein.radius_of_gyration())
    
    return np.array(rg_values)

# Calculate for all three systems
rg_wt = calculate_rg("wt_protein_only.gro", "wt_protein_only.xtc")
rg_chim = calculate_rg("chimera_protein_only.gro", "chimera_protein_only.xtc")
rg_af = calculate_rg("af_protein_only.gro", "Alpha_protein_only.xtc")

# Plotting with your specific color scheme
plt.figure(figsize=(10, 6))

# Define time axis based on the trajectory (500ns total)
def get_time(data):
    return np.linspace(0, 500, len(data))

plt.plot(get_time(rg_wt), rg_wt, label="Wild Type", color='blue', alpha=0.7)
plt.plot(get_time(rg_af), rg_af, label="AlphaFold Mutant", color='green', alpha=0.7)
plt.plot(get_time(rg_chim), rg_chim, label="Chimera Mutant", color='red', alpha=0.7)
plt.xlabel("Time (ns)")
plt.ylabel("Radius of Gyration (Å)")
plt.title("Protein Compactness ($R_g$) over 500ns")
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)
plt.show()

print(f"Average Rg - WT: {np.mean(rg_wt):.3f} Å")
print(f"Average Rg - AF: {np.mean(rg_af):.3f} Å")
print(f"Average Rg - CH: {np.mean(rg_chim):.3f} Å")
