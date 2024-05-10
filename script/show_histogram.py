import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams["font.size"] = 16

counts = []

with open("../build/braid_counts.txt", "r") as input_file:
    for _ in range(100):
        counts.append(int(input_file.readline().split()[1]))

fig, ax = plt.subplots(figsize = (7, 4))
plt.hist(counts, bins = range(1,101))
plt.xlim(1,101)
fig.subplots_adjust(left=0.08, right=0.95, bottom=0.16, top=0.95)
plt.xlabel("The number of homotopy classes")
plt.savefig("./histogram.pdf")
plt.show()
