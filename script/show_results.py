import matplotlib.pyplot as plt
import numpy as np
import matplotlib

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams["font.size"] = 16


PP_results=[]
PPvP_results=[]
CBS_results=[]
for i in range(100):
    with open("../build/CBS_opt_results2/"+"14-14-10-"+str(i)+".txt", "r") as res_file:
        res_file.readline()
        res_file.readline()
        res_file.readline()
        res_file.readline()
        cbs_c = float(res_file.readline())
        CBS_results.append(cbs_c)
    res = []
    with open("../build/PP_opt_results4/"+"14-14-10-"+str(i)+"-100.txt", "r") as res_file:
        res_file.readline()
        res_file.readline()
        res_file.readline()
        print(i)
        for t in range(100):
            res_file.readline()
            c = float(res_file.readline())
            res.append(c)
    PP_results.append(res)
    res = []
    with open("../build/PPvP_opt_results4/"+"14-14-10-"+str(i)+"-100.txt", "r") as res_file:
        res_file.readline()
        res_file.readline()
        res_file.readline()
        for t in range(100):
            res_file.readline()
            c = float(res_file.readline())
            res.append(c)
    PPvP_results.append(res)


fig, ax = plt.subplots(figsize = (7, 4))

mins = [[min(res[:i+1]) for i in range(len(res))] for res in PP_results]
mins=np.array(mins)
plt.plot(range(1,101), np.mean(mins, axis=0), label = "Ours")
#rates = [[1 if r<1.0 else 0 for r in rate] for rate in mins]
#rates = np.array(rates)
#plt.errorbar(range(100), np.mean(rates, axis=0), np.std(rates, axis=0), label = "PPHP")

mins = [[min(res[:i+1]) for i in range(len(res))] for res in PPvP_results]
mins=np.array(mins)
plt.plot(range(1,101), np.mean(mins, axis=0), label = "PPvP")
#rates = [[1 if r<1.0 else 0 for r in rate] for rate in mins]
#rates = np.array(rates)
#plt.errorbar(range(100), np.mean(rates, axis=0), np.std(rates, axis=0), label = "PPvP")

plt.plot(range(1,101), np.mean(cbs_c)*np.ones(100), label = "OO")

#rates = [[1 if r<1.0 else 0 for r in rate] for rate in mins]
#rates = np.array(rates)
#plt.errorbar(range(100), np.mean(rates, axis=0), np.std(rates, axis=0))
#plt.plot(range(100), np.mean(CBS_results)*np.ones(100))
#plt.plot(range(100), CBS_results)

plt.legend()
plt.xlim(1, 100)
plt.ylabel("Minimum cost")
plt.xlabel("The number of initial plans")
fig.subplots_adjust(left=0.14, right=0.95, bottom=0.16, top=0.95)
plt.savefig("./cost_graph.pdf")
plt.show()

