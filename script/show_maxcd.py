import matplotlib.pyplot as plt
import numpy as np
import matplotlib

def plot_region(ax, x, y, er, label, lw = 1, c = 0):
    color = matplotlib.cm.tab10(2+c)
    ax.plot(x, y, label = label, lw = lw, color = color)
    ax.fill_between(x, y-er, y+er, alpha = 0.3, color = color)

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams["font.size"] = 16

fig, ax = plt.subplots(figsize = (7, 4))

from matplotlib.ticker import ScalarFormatter

N = 10
K = 500

map_name = "empty-48-48"

max_cds1 = []
max_cds2 = []
for i in range(N):
    with open("../build/logs_Dyn_" + map_name + "/500-"+str(i)+"-100.txt", "r") as res_file:
        max_cd1 = []
        max_cd2 = []
        for _ in range(K):
            _, _, _, _, _, _, m1, m2 = map(float, res_file.readline().split())
            max_cd1.append(m1)
            max_cd2.append(m2)
        max_cds1.append(max_cd1)
        max_cds2.append(max_cd2)

max_cds1 = np.array(max_cds1)
max_cds2 = np.array(max_cds2)

plot_region(ax, range(2,K+1),np.log(max_cds1[:,1:]).mean(axis=0),np.log(max_cds1[:,1:]).std(axis=0), label = map_name, c = 0)

map_name = "den312d"

max_cds1 = []
max_cds2 = []
for i in range(N):
    with open("../build/logs_Dyn_" + map_name + "/500-"+str(i)+"-100.txt", "r") as res_file:
        max_cd1 = []
        max_cd2 = []
        for _ in range(K):
            _, _, _, _, _, _, m1, m2 = map(float, res_file.readline().split())
            max_cd1.append(m1)
            max_cd2.append(m2)
        max_cds1.append(max_cd1)
        max_cds2.append(max_cd2)

max_cds1 = np.array(max_cds1)
max_cds2 = np.array(max_cds2)

plot_region(ax, range(2,K+1),np.log(max_cds1[:,1:]).mean(axis=0),np.log(max_cds1[:,1:]).std(axis=0), label = map_name, c = 1)

map_name = "random-64-64-10"

max_cds1 = []
max_cds2 = []
for i in range(N):
    with open("../build/logs_Dyn_" + map_name + "/500-"+str(i)+"-100.txt", "r") as res_file:
        max_cd1 = []
        max_cd2 = []
        for _ in range(K):
            _, _, _, _, _, _, m1, m2 = map(float, res_file.readline().split())
            max_cd1.append(m1)
            max_cd2.append(m2)
        max_cds1.append(max_cd1)
        max_cds2.append(max_cd2)

max_cds1 = np.array(max_cds1)
max_cds2 = np.array(max_cds2)

plot_region(ax, range(2,K+1),np.log(max_cds1[:,1:]).mean(axis=0),np.log(max_cds1[:,1:]).std(axis=0), label = map_name, c = 2)

plt.xlim(1.9, K+1)
plt.ylabel("Log. max. coord.")
plt.xlabel("Number of agents")

plt.legend()

fig.subplots_adjust(left=0.10, right=0.97, bottom=0.16, top=0.95)
plt.savefig("./max_cd.pdf")
plt.show()
