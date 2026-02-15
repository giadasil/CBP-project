import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import f_oneway

file1 = pd.read_csv("Path\\to\\hbonds_chimera.dat", sep='\s+', names=["time", "num_hbonds"])
file2 = pd.read_csv("Path\\to\\hbonds_alpha.dat", sep='\s+', names=["time", "num_hbonds"])
file3 = pd.read_csv("Path\\to\\hbonds_wt.dat", sep='\s+', names=["time", "num_hbonds"])
print(file1.head())

#Check distribution normality
plt.subplot(1,3,1)
plt.hist(file1["num_hbonds"], bins=20, alpha=0.5, label="Chimera")
plt.subplot(1,3,2)
plt.hist(file2["num_hbonds"], bins=20, alpha=0.5, label="AlphaFold")
plt.subplot(1,3,3)
plt.hist(file3["num_hbonds"], bins=20, alpha=0.5, label="Wild Type")
plt.suptitle("Distribution of Hydrogen Bonds", size=20)
plt.show()

#plot boxplots of hbonds distribution
plt.subplot(1,3,1)
box3 = file3.boxplot(column="num_hbonds",rot=90, figsize=(3, 6))
plt.yticks(np.arange(6, 43, 0.5))
plt.title("Dist. of hBonds in Wild Type", size=20)
plt.ylabel("Number of Hydrogen Bonds", size=15)
plt.subplot(1,3,2)
box2 = file2.boxplot(column="num_hbonds",rot=90, figsize=(3, 6))
plt.yticks(np.arange(6, 43, 0.5))
plt.title("Dist. of hBonds in AlphaFold", size=20)
plt.ylabel("Number of Hydrogen Bonds", size=15)
plt.subplot(1,3,3)
box1 = file1.boxplot(column="num_hbonds",rot=90, figsize=(3, 6))
plt.title("Dist. of hBonds in Chimera", size=20)
plt.ylabel("Number of Hydrogen Bonds", size=15)
#plt.xlabel("Simulation Time (ns)", size=15)
plt.yticks(np.arange(6, 43, 0.5))
plt.show()

#Anova analysis
f_statistic, p_value = f_oneway(file1["num_hbonds"], file2["num_hbonds"])

print(f"F-statistic: {f_statistic}")
print(f"P-value: {p_value}")


