import matplotlib.pyplot as plt
import numpy as np
import matplotlib

def plot_region(ax, x, y, er, label, lw = 1):
    ax.plot(x, y, label = label, lw = lw)
    ax.fill_between(x, y-er, y+er, alpha = 0.3)

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams["font.size"] = 16

fig = plt.figure(figsize = (14, 4))

N = 10
K = 500

ax = fig.add_subplot(1,3,1)

runtimes=[]
closed_counts = []
all_counts = []

map_name = "empty-48-48"

for i in range(N):
    with open("../build/logs_Dyn_" + map_name + "/500-"+str(i)+"-100.txt", "r") as res_file:
        runtime = []
        length = []
        closed_count = []
        all_count = []
        _r = 0
        _c = 0
        for j in range(K):
            r, _, a, _, c, _, _, _ = map(int, res_file.readline().split())
            runtime.append(r)
            closed_count.append(c)
            all_count.append(a)
            _r = r
            _c = c
        runtimes.append(runtime)
        closed_counts.append(closed_count)
        all_counts.append(all_count)

runtimes = np.array(runtimes)/1000
plot_region(ax, range(2,K+1), runtimes[:,1:].mean(axis=0),runtimes[:,1:].std(axis=0)/np.sqrt(10), label = 'Dyn', lw =2)
closed_counts=np.array(closed_counts)
all_counts=np.array(all_counts)
print(np.log(runtimes[:,499].mean()/runtimes[:,449].mean())/np.log(5/4.5))

runtimes=[]
for i in range(N):
    with open("../build/logs_HR_" + map_name + "/500-"+str(i)+"-100.txt", "r") as res_file:
        runtime = []
        for _ in range(K):
            try:
                r, _, _, _, _, _, _, _ = map(int, res_file.readline().split())
            except:
                break
            if r >= 72000000:
                break
            runtime.append(r)
        runtimes.append(runtime)

_K = min([len(runtime) for runtime in runtimes])
runtimes = np.array([runtime[:_K] for runtime in runtimes])/1000
plot_region(ax, range(2,_K+1), runtimes[:,1:].mean(axis=0),runtimes[:,1:].std(axis=0)/np.sqrt(10), label = 'HR')
print(np.log(runtimes[:,49].mean()/runtimes[:,39].mean())/np.log(50/40))


plt.xlim(1.9, K+1)
plt.ylim(1.,24000)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylabel("Runtime (s)")
ax.set_xlabel("Number of agents")
ax.set_title("empty-48-48")
plt.legend()

map_name = "den312d"

ax = fig.add_subplot(1,3,2)

runtimes=[]
closed_counts = []
all_counts = []

for i in range(N):
    with open("../build/logs_Dyn_" + map_name + "/500-"+str(i)+"-100.txt", "r") as res_file:
        runtime = []
        length = []
        closed_count = []
        all_count = []
        _r = 0
        _c = 0
        for j in range(K):
            r, _, a, _, c, _, _, _ = map(int, res_file.readline().split())
            runtime.append(r)
            closed_count.append(c)
            all_count.append(a)
            _r = r
            _c = c
        runtimes.append(runtime)
        closed_counts.append(closed_count)
        all_counts.append(all_count)

runtimes = np.array(runtimes)/1000
plot_region(ax, range(2,K+1), runtimes[:,1:].mean(axis=0),runtimes[:,1:].std(axis=0)/np.sqrt(10), label = 'Dyn', lw =2)
closed_counts=np.array(closed_counts)
all_counts=np.array(all_counts)
print(np.log(runtimes[:,499].mean()/runtimes[:,449].mean())/np.log(5/4.5))

runtimes=[]
for i in range(N):
    with open("../build/logs_HR_" + map_name + "/500-"+str(i)+"-100.txt", "r") as res_file:
        runtime = []
        for _ in range(K):
            try:
                r, _, _, _, _, _, _, _ = map(int, res_file.readline().split())
            except:
                break
            if r >= 72000000:
                break
            runtime.append(r)
        runtimes.append(runtime)

_K = min([len(runtime) for runtime in runtimes])
runtimes = np.array([runtime[:_K] for runtime in runtimes])/1000
plot_region(ax, range(2,_K+1), runtimes[:,1:].mean(axis=0),runtimes[:,1:].std(axis=0)/np.sqrt(10), label = 'HR')
print(np.log(runtimes[:,49].mean()/runtimes[:,39].mean())/np.log(5/4))


plt.xlim(1.9, K+1)
plt.ylim(1.,24000)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylabel("Runtime (s)")
ax.set_xlabel("Number of agents")
ax.set_title("den312d")
plt.legend()


map_name = "random-64-64-10"

ax = fig.add_subplot(1,3,3)

runtimes=[]
closed_counts = []
all_counts = []

for i in range(N):
    with open("../build/logs_Dyn_" + map_name + "/500-"+str(i)+"-100.txt", "r") as res_file:
        runtime = []
        length = []
        closed_count = []
        all_count = []
        _r = 0
        _c = 0
        for j in range(K):
            r, _, a, _, c, _, _, _ = map(int, res_file.readline().split())
            runtime.append(r)
            closed_count.append(c)
            all_count.append(a)
            _r = r
            _c = c
        runtimes.append(runtime)
        closed_counts.append(closed_count)
        all_counts.append(all_count)

runtimes = np.array(runtimes)/1000
plot_region(ax, range(2,K+1), runtimes[:,1:].mean(axis=0),runtimes[:,1:].std(axis=0)/np.sqrt(10), label = 'Dyn', lw =2)
closed_counts=np.array(closed_counts)
all_counts=np.array(all_counts)
print(np.log(runtimes[:,499].mean()/runtimes[:,449].mean())/np.log(5/4.5))

runtimes=[]
for i in range(N):
    with open("../build/logs_HR_" + map_name + "/500-"+str(i)+"-100.txt", "r") as res_file:
        runtime = []
        for _ in range(K):
            try:
                r, _, _, _, _, _, _, _ = map(int, res_file.readline().split())
            except:
                break
            if r >= 72000000:
                break
            runtime.append(r)
        runtimes.append(runtime)

_K = min([len(runtime) for runtime in runtimes])
runtimes = np.array([runtime[:_K] for runtime in runtimes])/1000
plot_region(ax, range(2,_K+1), runtimes[:,1:].mean(axis=0),runtimes[:,1:].std(axis=0)/np.sqrt(10), label = 'HR')


plt.xlim(1.9, K+1)
plt.ylim(1.,24000)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylabel("Runtime (s)")
ax.set_xlabel("Number of agents")
ax.set_title("random-64-64-10")
plt.legend()

fig.subplots_adjust(left=0.06, right=0.99, bottom=0.16, top=0.90, wspace = 0.22)

plt.savefig("./runtime.pdf")
plt.show()
