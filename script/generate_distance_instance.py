#!/usr/bin/python3

import sys
import random
import itertools, copy

h = int(sys.argv[1])
w = int(sys.argv[2])
n = int(sys.argv[3])
seed = int(sys.argv[4])
random.seed(seed)

print(h,w,n)

states = []
pool = set()
for i in range(1,h-1):
    for j in range(1,w-1):
        pool.add((i,j))
for _ in range(2*n):
    s = random.sample(pool, 1)[0]
    for x in range(-1,2):
        for y in range(-1,2):
            pool.discard((s[0]+x, s[1]+y))
    states.append(s)
starts = states[::2]
goals = states[1::2]
for s,g in zip(starts, goals):
    print(s[0], s[1], g[0], g[1])
