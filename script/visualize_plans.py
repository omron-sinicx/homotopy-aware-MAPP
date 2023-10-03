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

time_span = float(sys.argv[2])
fig, ax = plt.subplots(figsize = (4,4))

#sin1 = math.sin(1.)
#cos1 = math.cos(1.)
sin1=0.
cos1=1.
with open(sys.argv[1], "r") as plan_file:
    width, height, number_of_agents, number_of_plans = map(int, plan_file.readline().split())
    ax.set_xlim(0, cos1*width+sin1*height)
    ax.set_ylim(-sin1*width, cos1*height)
    plans=[]
    agent_size = 0.5

    for _ in range(number_of_plans):
        line = plan_file.readline().split()
        makespan = int(line[0])
        if len(line)>1:
            delta = float(line[1])
        else:
            delta = 1.
        routes = []
    
        for _ in range(makespan+1):
            
            pos = list(map(float,plan_file.readline().split()))
            routes.append(pos)

        plans.append(routes)
        plan_file.readline()

    routes = plans[int(sys.argv[3])]
    agent_circles = []

    colors = plt.cm.get_cmap('jet', number_of_agents)
    
    def init():
        global steps
        global routes
            
        steps = number_of_agents * [0]
        init_patches = []
        
        for agent in range(number_of_agents):
            position = routes[-1][2*agent:2*agent+2]
            position += np.array([0.5, 0.5])

            init_patches.append(ax.add_patch(patches.Circle(position, radius = agent_size,fill = False, linewidth = 0.3, linestyle = '-', color = colors(agent), alpha = 0.5)))
        
        return []
    
    agent_circles = []
    
    def plot(frame):
        global routes
        global agent_circles
        for circle in agent_circles:
            circle.remove()

        t = time_span * frame / delta

        agent_circles = []
        step = math.floor(t)
        for agent in range(number_of_agents):
            if step >= len(routes)-1:
                position = np.array(routes[-1][2*agent:2*agent+2])
            else:
                r = t - math.floor(t)
                position = (1-r) * np.array(routes[step][2*agent:2*agent+2]) + r * np.array(routes[step+1][2*agent:2*agent+2])
            
            position += np.array([0.5, 0.5])
            agent_circles.append(ax.add_patch(patches.Circle([cos1*position[0]+sin1*position[1],-sin1*position[0]+cos1*position[1]], radius = agent_size, color = colors(agent))))
        return agent_circles
    
    ani = animation.FuncAnimation(fig, plot, frames = range(math.ceil(delta*(len(routes)-1) / time_span)+10), interval=100, init_func = init, blit = True)
    plt.show()
    if len(sys.argv) > 4:
        ani.save(sys.argv[4], writer='ffmpeg')
        #ani.save(sys.argv[4], writer='imagemagick')
