#!/usr/bin/env python3

import sys

with open(sys.argv[1], "r") as plan_file:
    width, height, number_of_agents, number_of_plans = map(int, plan_file.readline().split())
    
    plans=[]
    agent_size = 0.5

    for _ in range(number_of_plans):
        line = plan_file.readline().split()
        makespan = int(line[0])
        delta = 1.
        routes = []
    
        for _ in range(makespan+1):
            
            pos = list(map(int,plan_file.readline().split()))
            routes.append([(pos[2*i], pos[2*i+1]) for i in range(number_of_agents)])

        for i in range(number_of_agents):
            for j in range(i+1, number_of_agents):
                for t in range(makespan+1):
                    assert routes[t][i] != routes[t][j]
                for t in range(makespan):
                    #assert routes[t][i] != routes[t+1][j]
                    #assert routes[t+1][i] != routes[t][j]
                    assert routes[t][i] != routes[t+1][j] or routes[t+1][i] != routes[t][j]
        plans.append(routes)
        plan_file.readline()
