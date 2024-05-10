import matplotlib.pyplot as plt
import numpy as np
import matplotlib

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams["font.size"] = 16

fig, ax = plt.subplots(figsize = (5, 4))

for N in [3, 5, 10, 100, 1000]:
    file_path = "../build/cd_values/"+str(N)+'.txt'
    data = []
    with open(file_path, "r") as res_file:
        for _ in range(100):
            ms = []
            for _ in range(100):
                m = int(res_file.readline())
                ms.append(np.log(np.float128(m)))
            data.append(ms)
    data = np.array(data)
    plt.errorbar(range(10,1001,10), data.mean(axis=0), data.std(axis=0), label="n = "+str(N))

from matplotlib.ticker import ScalarFormatter

plt.legend()
plt.ylabel("Log. max. coord.")
plt.xlabel("l/n")
fig.subplots_adjust(left=0.16, right=0.95, bottom=0.16, top=0.95)
plt.savefig("./cd_values.pdf")
plt.show()
