import snap
import random
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


def runSIS (G, init_adopter, p, tl):
    number_node = G.GetNodes()
    number_edge = G.GetEdges()

    initial_list = random.sample(range(number_node), init_adopter)

    # data to return
    infectious_num = []
    removed_num = []
    final_states = []
    numberlist = [init_adopter]

    # 0: Susceptible, 1:Infectious
    labels = snap.TIntStrH()
    for NI in G.Nodes():
        nid = NI.GetId()
        state = 0
        if nid in initial_list:
            G.AddIntAttrDatN(nid, 1, "state")
            G.AddIntAttrDatN(nid, 1, "infected")
            state = 1
        else:
            G.AddIntAttrDatN(nid, 0, "state")
            G.AddIntAttrDatN(nid, 0, "infected")
            state = 0
        G.AddIntAttrDatN(nid, 0, "step")


    number_R = number_node - len(initial_list)
    number_echo = 0
    #number_I = 0
    #while number_R != number_node:
    while numberlist[-1] < number_node/2:
        print(numberlist[-1])
        for NI in G.Nodes():
            nid = NI.GetId()
            ival = G.GetIntAttrDatN(nid, "state")
            # if it's infectious, it will infect its neighbors.
            if ival == 1:
                for nid1 in NI.GetOutEdges():
                    att_value = G.GetIntAttrDatN(nid1, "state")

                    # with probability of p=0.2, a susceptible node will be infected.
                    if (np.random.choice(np.arange(0,2), p=[p, 1-p]) == 0) and (att_value == 0):
                        # The state 3 is temporary in order to avoid repeated infection.
                        G.AddIntAttrDatN(nid1, 3, "state")
                        G.AddIntAttrDatN(nid1, 1, "infected")

                # update the step of infection.
                step_this = G.GetIntAttrDatN(nid, "step")
                step_this += 1
                G.AddIntAttrDatN(nid, step_this, "step")
                if step_this == tl:
                    G.AddIntAttrDatN(nid, 4, "state")
                    G.AddIntAttrDatN(nid, 0, "step")

        # Count the number of susceptible or removed nodes.
        number_R_temp = 0
        for NI in G.Nodes():
            nid = NI.GetId()
            ival = G.GetIntAttrDatN(nid, "state")
            if ival == 3: # The nodes should be updated as infectious now.
                G.AddIntAttrDatN(nid, 1, "state")
                ival = 1
            elif ival == 4:
                G.AddIntAttrDatN(nid, 0, "state")
                ival = 0
            if ival == 0:
                number_R_temp += 1
        number_R = number_R_temp
        number_echo += 1
        infectious_num.append(number_node - number_R)
        #print(number_R)

        # Count the number of people that have been infected
        number_I = 0
        for NI in G.Nodes():
            nid = NI.GetId()
            infected = G.GetIntAttrDatN(nid, "infected")
            if infected == 1:
                number_I += 1
        numberlist.append(number_I)


    for NI in G.Nodes():
        nid = NI.GetId()
        ival = G.GetIntAttrDatN(nid, "state")
        step_number = G.GetIntAttrDatN(nid, "step")
        final_states.append((nid, ival, step_number))
        # print(nid, ival, step_number)

    print("One round of simulation completes.")
    return infectious_num, final_states



# ---------RUN MODELS & COLLECT DATA----------
if len(sys.argv) != 2:
    print("Usage: <figure #>\n1. tl - finish time\n2. p - finish time\n3. # nodes in each state - step number\n")

G = snap.LoadEdgeList(snap.PNEANet, "soc-Epinions1.txt", 0, 1, '\t')
number_node = G.GetNodes()

init_adopter = 10
p = 0.2
tl = 2

# tl - finish time fig
if sys.argv[1] == '1' or sys.argv[1] == 'all':
    finish_steps = []
    for tl in range(1, 6):
        infectious_num, final_states = runSIS(G, init_adopter, p, tl)
        finish_steps.append(len(infectious_num))
    remove_interval = [x for x in range(1, 6)]
    print remove_interval, finish_steps
    fig, ax = plt.subplots()
    ax.plot(remove_interval, finish_steps, '.-')

    ax.set(xlabel='time interval to enter remove state', ylabel='step that the infection stops',
           title='SIS - time interval and stop time')
    plt.xticks(remove_interval)
    plt.yticks(finish_steps)
    ax.grid()

    fig.savefig("SIS-tl-finish_step.png")

# p - finish time fig
if sys.argv[1] == '2' or sys.argv[1] == 'all':
    finish_steps = []
    for i in range(1, 10):
        infectious_num, final_states = runSIS(G, init_adopter, float(i)/10, tl)
        finish_steps.append(len(infectious_num))
    contagion_prob  = [float(x)/10 for x in range(1, 10)]
    print contagion_prob, finish_steps
    fig, ax = plt.subplots()
    ax.plot(contagion_prob, finish_steps, '.-')

    ax.set(xlabel='Contagion Probability', ylabel='step that the infection stops',
           title='SIS - contagion probability and stop time')
    plt.xticks(contagion_prob)
    plt.yticks(finish_steps)
    ax.grid()

    fig.savefig("SIS-p-finish_step.png")

# infected/susceptible - step number fig
if sys.argv[1] == '3' or sys.argv[1] == 'all':
    infectious_num, final_states = runSIS(G, init_adopter, p, tl)
    print(infectious_num)
    steps = [x for x in range (0, len(infectious_num)+1)]
    zero = [0]
    infectious_num = np.array(zero + infectious_num)
    susceptible_num = np.array([(number_node - infectious_num[x]) for x in steps])

    fig, ax = plt.subplots()
    width = 0.2
    p1 = ax.bar(steps, susceptible_num, width, color='C0')
    p2 = ax.bar(steps, infectious_num, width, bottom=susceptible_num, color='C1')

    ax.set(xlabel='step number', ylabel='# nodes',
           title='SIS - evolution of the number of nodes in each state')
    plt.xticks(steps)
    ax.legend((p1[0], p2[0]), ('susceptible nodes', 'infectious nodes'))

    fig.savefig("SIS-step-#state.png")