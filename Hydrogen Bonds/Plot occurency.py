import matplotlib.pyplot as plt
import pandas as pd

#Select the file hbonds_details.dat to read the data from
file = pd.read_csv("Path\\to\\hbonds_details.dat",skiprows=1, sep='\s+')


#Occupancy is given as a string with percentage sign
#Convert occupancy to a float
print(file.head())
file["occupancy"] = (
    file["occupancy"]
    .str.replace("%", "", regex=False)
    .astype(float)
)

#Sort to have most present H-bonds at the top of the list
file_sorted = file.sort_values("occupancy", ascending=False)

#Select the top 10 most present H-bonds and plot them in a horizontal bar chart
fig, ax = plt.subplots()
hbar = ax.barh(
    file_sorted["donor"].head(10) + " – " + file_sorted["acceptor"].head(10),
    file_sorted["occupancy"].head(10),
    
)
labels = x = range(0, 130, 10)
plt.xticks(x, labels)
ax.tick_params(axis="y", labelsize=12) 
plt.title("Top 10 Hydrogen Bonds in Wild Type", size=20)
ax.set_xlabel("H-bond occupancy (%)", size=15)
ax.bar_label(hbar, fmt='%.2f', size=15)
plt.show()