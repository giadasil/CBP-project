import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

u_wt = mda.Universe("wt_protein_only.gro", "wt_protein_only.xtc")
u_af = mda.Universe("af_protein_only.gro", "Alpha_protein_only.xtc")
u_ch = mda.Universe("chimera_protein_only.gro", "chimera_protein_only.xtc")

# 1. Calculate start frame 
total_frames = len(u_wt.trajectory)
start_frame = total_frames - 9000

def get_equilibrated_rg(u, start):
    protein = u.select_atoms("name CA")
    rg_values = []
    # Loop only through the last 9000 frames
    for ts in u.trajectory[start:]:
        rg_values.append(u.atoms.radius_of_gyration())
    return np.array(rg_values)

# 2. Run analysis
rg_wt_equil = get_equilibrated_rg(u_wt, start_frame)
rg_af_equil = get_equilibrated_rg(u_af, start_frame)
rg_ch_equil = get_equilibrated_rg(u_ch, start_frame)

# 3. Print final comparison table
print(f"{'System':<20} | {'Mean Rg (Å)':<15} | {'Std Dev (Å)':<10}")
print("-" * 50)
systems_data = [
    ("Wild Type", rg_wt_equil),
    ("AlphaFold Mutant", rg_af_equil),
    ("Chimera Mutant", rg_ch_equil)
]

for name, data in systems_data:
    print(f"{name:<20} | {np.mean(data):.3f} Å      | {np.std(data):.3f}")

plt.figure(figsize=(10, 5))
time_equil = np.array([ts.time for ts in u_wt.trajectory[start_frame:]])
time_ns = time_equil / 1000.0

plt.plot(time_ns, rg_wt_equil, label="WT", color='blue', alpha=0.6)
plt.plot(time_ns, rg_af_equil, label="AF Mutant", color='green', alpha=0.6)
plt.plot(time_ns, rg_ch_equil, label="Chimera Mutant", color='red', alpha=0.6)

plt.xlabel("Time (ns)")
plt.ylabel("Radius of Gyration (Å)")
plt.title("Protein Compactness (Last 9000 Frames)")
plt.legend()
plt.show()
