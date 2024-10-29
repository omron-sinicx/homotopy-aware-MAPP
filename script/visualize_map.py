#!/usr/bin/env python3

# Copyright (c) 2022 OMRON SINIC X Corporation
# Author: Kazumi Kasaura


import matplotlib.pyplot as plt 
import matplotlib.patches as patches
import sys
import matplotlib.animation as animation
import math
import numpy as np
import matplotlib.path

obstacles = []
fig, ax = plt.subplots(figsize = (8,8))

with open(sys.argv[1], "r") as map_file:
    map_file.readline()
    height = int(map_file.readline().split()[1])
    width = int(map_file.readline().split()[1])
    map_file.readline()
    ax.set_xlim(0, width)
    ax.set_ylim(0, height)
    # #rows lines with the map
    my_map = [[False for _ in range(height)] for _ in range(width)]
    for r in range(height):
        line = map_file.readline()
        for c in range(width):
            if line[c] == '.' or line[c] == 'G':
                my_map[c][height - r - 1] = False
            else:
                obstacles.append((c, height - r - 1))
                my_map[c][height - r - 1] = True

for x,y in obstacles:
    ax.add_patch(patches.Rectangle((x,y), 1., 1., color = 'black'))

ax.set_xticklabels([])
ax.set_yticklabels([])

ax.tick_params(left = False, bottom = False)

fig.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.99)

plt.savefig(sys.argv[2]+".pdf")

plt.show()
