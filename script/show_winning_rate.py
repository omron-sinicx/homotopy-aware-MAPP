import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import sys

def plot_region(ax, x, y, er, label, lw = 1):
    ax.plot(x, y, label = label, lw = lw)
    ax.fill_between(x, y-er, y+er, alpha = 0.3)

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams["font.size"] = 16

fig = plt.figure(figsize = (14, 4))

map_name = "empty-14-14"

PP_results=[]
PPvP_results=[]
CBS_results=[]

PP_wins = []
PPvP_wins = []
CBS_wins = []

for i in range(27):
    with open("../build/" + map_name + "_ICBS_opt_results/"+"10-"+str(i)+".txt", "r") as res_file:
        res_file.readline()
        res_file.readline()
        res_file.readline()
        res_file.readline()
        cbs_c = float(res_file.readline())
        CBS_results.append(cbs_c)
    PP_res = []
    with open("../build/" + map_name + "_PP_opt_results/"+"10-"+str(i)+"-100.txt", "r") as res_file:
        res_file.readline()
        res_file.readline()
        res_file.readline()
        for t in range(100):
            res_file.readline()
            c = float(res_file.readline())
            PP_res.append(c)
    PP_results.append(PP_res)
    PPvP_res = []
    with open("../build/" + map_name + "_PPvP_opt_results/"+"10-"+str(i)+"-100.txt", "r") as res_file:
        res_file.readline()
        res_file.readline()
        res_file.readline()
        for t in range(100):
            res_file.readline()
            c = float(res_file.readline())
            PPvP_res.append(c)
    PPvP_results.append(PPvP_res)

    cbs_win = np.zeros(100)
    pp_win = np.zeros(100)
    ppvp_win = np.zeros(100)
    for t in range(100):
        pp_c = min(PP_res[:t+1])
        ppvp_c = min(PPvP_res[:t+1])
        pp_win[t] = 1.0 if pp_c <= ppvp_c and pp_c <= cbs_c else 0.0
        ppvp_win[t] = 1.0 if ppvp_c <= pp_c and ppvp_c <= cbs_c else 0.0
        cbs_win[t] = 1.0 if cbs_c <= pp_c and cbs_c <= ppvp_c else 0.0

    PP_wins.append(pp_win)
    PPvP_wins.append(ppvp_win)
    CBS_wins.append(cbs_win)

ax = fig.add_subplot(1,2,1)

PP_wins = np.array(PP_wins)
PPvP_wins = np.array(PPvP_wins)
CBS_wins = np.array(CBS_wins)

plot_region(ax, range(1, 101), np.mean(PP_wins, axis=0), np.std(PP_wins, axis=0)/10, label = "Ours", lw = 2)
plot_region(ax, range(1, 101), np.mean(PPvP_wins, axis=0), np.std(PPvP_wins, axis=0)/10, label = "PPvP")
plot_region(ax, range(1, 101), np.mean(CBS_wins, axis=0), np.std(CBS_wins, axis=0)/10, label = "OO")

plt.legend()
ax.set_xlim(1, 100)
ax.set_ylim(0., 0.7)
ax.set_title("Empty")
ax.set_ylabel("Winning rate")
ax.set_xlabel("Number of initial plans")


map_name = "obstacle-14-14"

PP_results=[]
PPvP_results=[]
CBS_results=[]

PP_wins = []
PPvP_wins = []
CBS_wins = []

for i in range(27):
    with open("../build/" + map_name + "_ICBS_opt_results/"+"10-"+str(i)+".txt", "r") as res_file:
        res_file.readline()
        res_file.readline()
        res_file.readline()
        res_file.readline()
        cbs_c = float(res_file.readline())
        CBS_results.append(cbs_c)
    PP_res = []
    with open("../build/" + map_name + "_PP_opt_results/"+"10-"+str(i)+"-100.txt", "r") as res_file:
        res_file.readline()
        res_file.readline()
        res_file.readline()
        for t in range(100):
            res_file.readline()
            c = float(res_file.readline())
            PP_res.append(c)
    PP_results.append(PP_res)
    PPvP_res = []
    with open("../build/" + map_name + "_PPvP_opt_results/"+"10-"+str(i)+"-100.txt", "r") as res_file:
        res_file.readline()
        res_file.readline()
        res_file.readline()
        for t in range(100):
            res_file.readline()
            c = float(res_file.readline())
            PPvP_res.append(c)
    PPvP_results.append(PPvP_res)

    cbs_win = np.zeros(100)
    pp_win = np.zeros(100)
    ppvp_win = np.zeros(100)
    for t in range(100):
        pp_c = min(PP_res[:t+1])
        ppvp_c = min(PPvP_res[:t+1])
        pp_win[t] = 1.0 if pp_c <= ppvp_c and pp_c <= cbs_c else 0.0
        ppvp_win[t] = 1.0 if ppvp_c <= pp_c and ppvp_c <= cbs_c else 0.0
        cbs_win[t] = 1.0 if cbs_c <= pp_c and cbs_c <= ppvp_c else 0.0

    PP_wins.append(pp_win)
    PPvP_wins.append(ppvp_win)
    CBS_wins.append(cbs_win)

ax = fig.add_subplot(1,2,2)

PP_wins = np.array(PP_wins)
PPvP_wins = np.array(PPvP_wins)
CBS_wins = np.array(CBS_wins)

plot_region(ax, range(1, 101), np.mean(PP_wins, axis=0), np.std(PP_wins, axis=0)/10, label = "Ours", lw = 2)
plot_region(ax, range(1, 101), np.mean(PPvP_wins, axis=0), np.std(PPvP_wins, axis=0)/10, label = "PPvP")
plot_region(ax, range(1, 101), np.mean(CBS_wins, axis=0), np.std(CBS_wins, axis=0)/10, label = "OO")

plt.legend()
ax.set_xlim(1, 100)
ax.set_ylim(0., 0.7)
ax.set_title("Obstacle")
ax.set_ylabel("Winning rate")
ax.set_xlabel("Number of initial plans")



fig.subplots_adjust(left=0.07, right=0.98, bottom=0.16, top=0.90, wspace = 0.19)

plt.savefig("./wr_graph.pdf")
plt.show()

