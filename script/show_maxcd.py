import matplotlib.pyplot as plt
import numpy as np
import matplotlib

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams["font.size"] = 16

N = 10
K = 800

max_cds1 = []
max_cds2 = []
for i in range(N):
    with open("../build/logs7/"+"122-122-800-"+str(i)+"-100.txt", "r") as res_file:
        max_cd1 = []
        max_cd2 = []
        for _ in range(K):
            _, _, _, _, _, m1, m2 = map(int, res_file.readline().split())
            max_cd1.append(m1)
            max_cd2.append(m2)
        max_cds1.append(max_cd1)
        max_cds2.append(max_cd2)

max_cds1 = np.array(max_cds1)
max_cds2 = np.array(max_cds2)

fig, ax = plt.subplots(figsize = (5, 4))

from matplotlib.ticker import ScalarFormatter
plt.xlim(1.9, K+1)
plt.ylabel("Log. max. coord.")
plt.xlabel("The number of agents")

plt.errorbar(range(2,K+1),np.log(max_cds1[:,1:]).mean(axis=0),np.log(max_cds1[:,1:]).std(axis=0))

fig.subplots_adjust(left=0.14, right=0.96, bottom=0.16, top=0.95)
plt.savefig("./max_cd.pdf")
plt.show()
