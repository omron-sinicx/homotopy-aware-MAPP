import matplotlib.pyplot as plt
import numpy as np
import matplotlib

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams["font.size"] = 16


runtimes=[]
lengths=[]
for i in range(10):
    with open("../build/logs/"+"29-29-40-"+str(i)+"-100.txt", "r") as res_file:
        runtime = []
        length = []
        for _ in range(40):
            r, l = map(int, res_file.readline().split())
            runtime.append(r)
            length.append(l)
        runtimes.append(runtime)
        lengths.append(length)

runtimes = np.array(runtimes)/1000.
lengths = np.array(lengths)

fig, ax = plt.subplots(figsize = (5, 4))


plt.xlim(1.9, 41)
plt.ylim(bottom=0.101, top=700)
ax.set_xscale('log')
ax.set_yscale('log')
#ax.set_aspect('equal')
plt.ylabel("Braid word length")
plt.xlabel("Number of agents")
plt.errorbar(range(2,41), lengths[:,1:].mean(axis=0),lengths[:,1:].std(axis=0))
fig.subplots_adjust(left=0.16, right=0.98, bottom=0.16, top=0.95)
plt.savefig("./length.pdf")
plt.show()

