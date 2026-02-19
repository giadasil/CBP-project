import MDAnalysis as mda
from MDAnalysis.analysis import pca, align
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
from sklearn.cluster import DBSCAN
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors

# 1. Load Universes
u_wt = mda.Universe("wt_protein_only.gro", "wt_protein_only.xtc")
u_af = mda.Universe("af_protein_only.gro", "Alpha_protein_only.xtc")
u_ch = mda.Universe("chimera_protein_only.gro", "chimera_protein_only.xtc")

# 2. Pre-align trajectories to the WT reference structure
ref_atoms = u_wt.select_atoms('name CA')
for u in [u_wt, u_af, u_ch]:
    alignment = align.AlignTraj(u, u_wt, select='name CA', in_memory=True).run()

# 3. Run PCA on the WT to define the "Reference Space"
pc_wt = pca.PCA(u_wt, select='name CA').run()

# Project all systems onto the first 3 WT axes
wt_3d = pc_wt.transform(u_wt.select_atoms('name CA'), n_components=3)
af_3d = pc_wt.transform(u_af.select_atoms('name CA'), n_components=3)
ch_3d = pc_wt.transform(u_ch.select_atoms('name CA'), n_components=3)

# 4. Data Preparation for the Last 900 Frames
def get_pca_df(projections, label, last_frames=900):
    df = pd.DataFrame(projections, columns=['PC1', 'PC2', 'PC3'])
    df['System'] = label
    return df.tail(last_frames)

df_wt = get_pca_df(wt_3d, 'Wild Type')
df_af = get_pca_df(af_3d, 'Gln104Ala AF')
df_ch = get_pca_df(ch_3d, 'Gln104Ala Chimera')

# Combine systems into a single dataset for clustering comparison
df_all = pd.concat([df_wt, df_af, df_ch]).reset_index(drop=True)

# 5. Scaling (Critical for DBSCAN)
scaler = StandardScaler()
X_scaled = scaler.fit_transform(df_all[['PC1', 'PC2', 'PC3']])

# 6. DBSCAN Clustering with Initialization
# Initialize score to 0 to avoid NameError
score = 0.0

# Try different eps (e.g., 0.3 instead of 0.4)
# Or adjust min_samples to be more or less restrictive
dbscan = DBSCAN(eps=0.3, min_samples=30)
df_all['Cluster'] = dbscan.fit_predict(X_scaled)

# 7. Silhouette Score Validation with Error Handling
mask = df_all['Cluster'] != -1
labels = df_all.loc[mask, 'Cluster']
n_clusters = len(set(labels))

if n_clusters > 1:
    score = silhouette_score(X_scaled[mask], labels)
    print(f"DBSCAN Results:")
    print(f"Clusters Detected: {n_clusters}")
    print(f"Noise Points:      {list(df_all['Cluster']).count(-1)}")
    print(f"Silhouette Score:  {score:.4f}")
else:
    print(f"Insufficient clusters found ({n_clusters}). Only one basin or all noise detected.")
    print("Try reducing 'eps' if you have only 1 cluster, or increasing it if you have too much noise.")

# 8. Final Interactive Plot
fig = px.scatter_3d(
    df_all, 
    x='PC1', y='PC2', z='PC3',
    color='Cluster',
    symbol='System',
    title=f"DBSCAN Clustering (Detected: {n_clusters} clusters, Score: {score:.3f})",
    opacity=0.7,
    color_continuous_scale='Viridis'
)
fig.update_traces(marker=dict(size=3))
fig.show()
