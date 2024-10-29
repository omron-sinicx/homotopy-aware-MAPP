import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams["font.size"] = 16

fig = plt.figure(figsize = (14, 4))

map_name = "empty-14-14"
counts = []

with open("../build/count_" + map_name + ".txt", "r") as input_file:
    for _ in range(100):
        counts.append(int(input_file.readline().split()[1]))

ax = fig.add_subplot(1,2,1)
ax.hist(counts, bins = range(1,101))
ax.set_xlim(1,101)
ax.set_ylim(0,20)
ax.set_title("Empty")
ax.set_xlabel("Number of homotopy classes")


map_name = "obstacle-14-14"
counts = []

with open("../build/count_" + map_name + ".txt", "r") as input_file:
    for _ in range(100):
        counts.append(int(input_file.readline().split()[1]))

ax = fig.add_subplot(1,2,2)
ax.hist(counts, bins = range(1,101))
ax.set_xlim(1,101)
ax.set_ylim(0,20)
ax.set_title("Obstacle")
ax.set_xlabel("Number of homotopy classes")

fig.subplots_adjust(left=0.04, right=0.98, bottom=0.16, top=0.90, wspace = 0.1)

plt.savefig("./histogram.pdf")
plt.show()
