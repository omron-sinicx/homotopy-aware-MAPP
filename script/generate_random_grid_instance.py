#!/usr/bin/python3

import sys
import random
import itertools, copy

h = int(sys.argv[1])
w = int(sys.argv[2])
n = int(sys.argv[3])
seed = int(sys.argv[4])
random.seed(seed)

print(h,w,n,-1)

states = list(itertools.product(range(h), range(w)))
random.shuffle(states)
starts = copy.deepcopy(states[:n])
random.shuffle(states)
goals = copy.deepcopy(states[:n])
for s,g in zip(starts, goals):
    print(s[0], s[1], g[0], g[1])
