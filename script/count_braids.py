import subprocess
import numpy as np
import sys
counts = []
map_name = sys.argv[1]
nob = sys.argv[2]

for i in range(100):
    ret = subprocess.run("../build/count_braids " + nob + " < ../build/" + map_name + "_PPvP_results/10-"+str(i)+"-100.txt", capture_output = True, text = True, shell = True)
    counts.append(int(ret.stdout))
    print(i, counts[-1])

print(np.mean(counts))
