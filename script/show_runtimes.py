import matplotlib.pyplot as plt
import numpy as np
import matplotlib

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams["font.size"] = 16

fig, ax = plt.subplots(figsize = (5, 4))

N = 10
K = 800

runtimes=[]
closed_counts = []
all_counts = []

for i in range(N):
    with open("../build/logs7/"+"122-122-800-"+str(i)+"-100.txt", "r") as res_file:
        runtime = []
        length = []
        closed_count = []
        all_count = []
        _r = 0
        _c = 0
        for _ in range(K):
            r, _, a, _, c, _, _ = map(int, res_file.readline().split())
            runtime.append(r)
            closed_count.append(c)
            all_count.append(a)
            _r = r
            _c = c
        runtimes.append(runtime)
        closed_counts.append(closed_count)
        all_counts.append(all_count)

runtimes = np.array(runtimes)/1000
plt.errorbar(range(2,K+1), runtimes[:,1:].mean(axis=0),runtimes[:,1:].std(axis=0), label = 'Dyn')
closed_counts=np.array(closed_counts)
all_counts=np.array(all_counts)

runtimes=[]
for i in range(N):
    with open("../build/logs8/"+"122-122-800-"+str(i)+"-100.txt", "r") as res_file:
        runtime = []
        for _ in range(K):
            try:
                r, _, _, _, _, _, _ = map(int, res_file.readline().split())
            except:
                break
            runtime.append(r)
        runtimes.append(runtime)

_K = min([len(runtime) for runtime in runtimes])
runtimes = np.array([runtime[:_K] for runtime in runtimes])/1000
plt.errorbar(range(2,_K+1), runtimes[:,1:].mean(axis=0),runtimes[:,1:].std(axis=0), label = 'HR')

plt.xlim(1.9, K+1)
ax.set_xscale('log')
ax.set_yscale('log')
plt.ylabel("Runtime (s)")
plt.xlabel("The number of agents")
plt.legend()
fig.subplots_adjust(left=0.16, right=0.98, bottom=0.16, top=0.95)
plt.savefig("./runtime.pdf")
plt.show()
