import matplotlib.pyplot as plt
import numpy as np
import matplotlib

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams["font.size"] = 16


runtimes=[]
lengths=[]
makespans = []
for i in range(10):
    with open("../build/logs/"+"29-29-40-"+str(i)+"-100.txt", "r") as res_file:
        runtime = []
        length = []
        makespan = []
        for _ in range(40):
            r, l, m = map(int, res_file.readline().split())
            runtime.append(r)
            length.append(l)
            makespan.append(m)
        runtimes.append(runtime)
        lengths.append(length)
        makespans.append(makespan)

runtimes = np.array(runtimes)/1000.
lengths = np.array(lengths)
makespans = np.array(makespans)

fig, ax = plt.subplots(figsize = (5, 4))

from matplotlib.ticker import ScalarFormatter
#plt.legend()
plt.xlim(1.9, 41)
#plt.ylim(16,43)
ax.set_xscale('log')
ax.set_yscale('log')
ax.yaxis.set_major_formatter(ScalarFormatter())
ax.yaxis.set_minor_formatter(ScalarFormatter())
#ax.set_aspect('equal')
plt.ylabel("Makespan")
plt.xlabel("Number of agents")
plt.errorbar(range(2,41),makespans[:,1:].mean(axis=0),makespans[:,1:].std(axis=0))
fig.subplots_adjust(left=0.16, right=0.98, bottom=0.16, top=0.95)
plt.savefig("./makespan.pdf")
plt.show()
