#!/usr/bin/python3

import sys
import random
import itertools, copy
import numpy as np
from queue import Queue

map_file_name = sys.argv[1]

n = int(sys.argv[2])

obstacles = []
with open(map_file_name, "r") as map_file:
    map_file.readline()
    h = int(map_file.readline().split()[1])
    w = int(map_file.readline().split()[1])
    map_file.readline()
    for y in range(h):
        s = map_file.readline()
        for x in range(w):
            if s[x] != '.' and s[x] != 'G':
                obstacles.append((x,h - y - 1))

seed = int(sys.argv[3])
random.seed(seed)

print(map_file_name)
print(n)

pool = set()
for i in range(w):
    for j in range(h):
        pool.add((i,j))
for obs in obstacles:
    pool.discard(obs)

starts= random.sample(pool, n)
goals = random.sample(pool, n)
for s,g in zip(starts, goals):
    print(s[0], s[1], g[0], g[1])
