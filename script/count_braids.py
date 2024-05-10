import subprocess
import numpy as np
counts = []
for i in range(100):
    ret = subprocess.run("../build/count_braids < ../build/PPvP_results3/14-14-10-"+str(i)+"-100.txt", capture_output = True, text = True, shell = True)
    counts.append(int(ret.stdout))
    print(i, counts[-1])

print(np.mean(counts))
