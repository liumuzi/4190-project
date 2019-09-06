import snap
import random
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


def runSIR (G, init_adopter, p, tl, initial_list=None):
    number_node = G.GetNodes()
    number_edge = G.GetEdges()
    if initial_list is None:
        initial_list = random.sample(range(number_node), init_adopter)

    # data to return
    infectious_num = []
    removed_num = []
    final_states = []

    # 0: Susceptible, 1:Infectious, 2: Removed
    labels = snap.TIntStrH()
    for NI in G.Nodes():
        nid = NI.GetId()
        state = 0
        if nid in initial_list:
            G.AddIntAttrDatN(nid, 1, "state")
            state = 1
        else:
            G.AddIntAttrDatN(nid, 0, "state")
            state = 0
        G.AddIntAttrDatN(nid, 0, "step")


    number_R = number_node - len(initial_list)
    number_echo = 0
    while number_R != number_node:
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

                # update the step of infection.
                step_this = G.GetIntAttrDatN(nid, "step")
                step_this += 1
                G.AddIntAttrDatN(nid, step_this, "step")
                if step_this == tl:
                    G.AddIntAttrDatN(nid, 2, "state")

        # Count the number of susceptible or removed nodes.
        number_R_temp = 0
        removed = 0
        for NI in G.Nodes():
            nid = NI.GetId()
            ival = G.GetIntAttrDatN(nid, "state")
            if ival == 3: # The nodes should be updated as infectious now.
                G.AddIntAttrDatN(nid, 1, "state")
                ival = 1

            if ival == 0 or ival == 2:
                number_R_temp += 1
            if ival == 2:
                removed += 1
        number_R = number_R_temp
        number_echo += 1
        infectious_num.append(number_node - number_R)
        removed_num.append(removed)
        #print(number_R)


    for NI in G.Nodes():
        nid = NI.GetId()
        ival = G.GetIntAttrDatN(nid, "state")
        step_number = G.GetIntAttrDatN(nid, "step")
        final_states.append((nid, ival, step_number))
        # print(nid, ival, step_number)

    print("One round of simulation completes.")
    return infectious_num, removed_num, final_states



# ---------RUN MODELS & COLLECT DATA----------
if len(sys.argv) != 2:
    print("Usage: <figure #>\n1. tl - finish time & removed #\n2. p - finish time & removed #\n3. # nodes in each state - step number\n4. initial adopter size - finish time & removed #\n5. largest & smallest initial adopters\n")

G = snap.LoadEdgeList(snap.PNEANet, "soc-Epinions1.txt", 0, 1, '\t')
number_node = G.GetNodes()

init_adopter = 10
p = 0.2
tl = 2

# tl - finish time fig
if sys.argv[1] == '1' or sys.argv[1] == 'all':
    finish_steps = [14,15,24,28,29]
    fnl_removed = [17113,25548,30851,34793,37606]
    # for tl in range(1, 6):
    #     infectious_num, removed_num, final_states = runSIR(G, init_adopter, p, tl)
    #     finish_steps.append(len(infectious_num))
    #     fnl_removed.append(removed_num[-1])
    remove_interval = [x for x in range(1, 6)]
    print remove_interval, finish_steps, fnl_removed
    fig, ax1 = plt.subplots()
    ax1.bar(remove_interval, fnl_removed, 0.4,color='y')
    ax1.set_ylabel('removed nodes', color='y')
    ax1.set_xlabel('Time Interval of infectious state')
    ax1.tick_params('y', colors='y')

    ax2 = ax1.twinx()
    ax2.plot(remove_interval, finish_steps, '.-')
    ax2.set_ylabel('step that the infection stops', color='C0')
    ax2.tick_params('y', colors='C0')
    ax2.set_title('SIR - time interval and stop time')
    plt.xticks(remove_interval)
    plt.yticks(finish_steps)
    fig.savefig("SIR-tl-finish_step.png")

# p - finish time & removed # fig
if sys.argv[1] == '2' or sys.argv[1] == 'all':
    finish_steps = []
    fnl_removed = []
    for i in range(1, 10):
        infectious_num, removed_num, final_states = runSIR(G, init_adopter, float(i)/10, tl)
        finish_steps.append(len(infectious_num))
        fnl_removed.append(removed_num[-1])
    contagion_prob  = [float(x)/10 for x in range(1, 10)]
    print contagion_prob, finish_steps, fnl_removed
    fig, ax1 = plt.subplots()
    ax1.bar(contagion_prob, fnl_removed, 0.05,color='y')
    ax1.set_ylabel('removed nodes', color='y')
    ax1.set_xlabel('Contagion Probability')
    ax1.tick_params('y', colors='y')

    ax2 = ax1.twinx()
    ax2.plot(contagion_prob, finish_steps, '.-')
    ax2.set_ylabel('step that the infection stops', color='C0')
    ax2.tick_params('y', colors='C0')
    ax2.set_title('SIR - contagion probability and stop time')
    plt.xticks(contagion_prob)
    plt.yticks(finish_steps)

    fig.savefig("SIR-p-finish_step.png")

# infected/removed/susceptible - step number fig
if sys.argv[1] == '3' or sys.argv[1] == 'all':
    infectious_num, removed_num, final_states = runSIR(G, init_adopter, p, tl)
    steps = [x for x in range (0, len(infectious_num)+1)]
    zero = [0]
    infectious_num = np.array(zero + infectious_num)
    removed_num = np.array(zero + removed_num)
    susceptible_num = np.array([(number_node - infectious_num[x] - removed_num[x]) for x in steps])

    fig, ax = plt.subplots()
    width = 0.2
    p1 = ax.bar(steps, susceptible_num, width)
    p2 = ax.bar(steps, infectious_num, width, bottom=susceptible_num)
    p3 = ax.bar(steps, removed_num, width, bottom=susceptible_num+infectious_num)

    ax.set(xlabel='step number', ylabel='# nodes',
           title='SIR - evolution of the number of nodes in each state')
    plt.xticks(steps)
    ax.legend((p1[0], p2[0], p3[0]), ('susceptible nodes', 'infectious nodes', 'removed nodes'))

    fig.savefig("SIR-step-#state.png")

# initial adopter size - finish time & removed #
if sys.argv[1] == '4' or sys.argv[1] == 'all':
    finish_steps = [19, 16, 16, 17, 15, 18, 15, 14, 14, 13, 14, 13, 15, 15, 13, 13, 13, 13, 13, 14, 14, 12, 13, 12, 15, 15, 12, 12, 12, 11, 14, 11, 11, 11, 13, 14, 12, 12, 11, 12, 13, 11, 13, 11, 11, 11, 11, 12, 11, 10]
    fnl_removed = [25280, 25699, 26019, 26331, 26642, 26936, 27203, 27426, 27846, 28381, 28666, 29061, 29175, 29392, 29801, 30237, 30434, 30896, 31032, 31278, 31713, 31811, 32159, 32659, 32695, 33188, 33438, 33652, 33811, 34551, 34715, 34756, 35156, 35467, 35681, 36218, 36333, 36642, 36937, 37324, 37457, 37703, 38013, 38327, 38543, 38989, 39255, 39645, 39974, 40150]
    # for init_adopter in range(10, 20000)[::400]:
    #     infectious_num, removed_num, final_states = runSIR(G, init_adopter, p, tl)
    #     finish_steps.append(len(infectious_num))
    #     fnl_removed.append(removed_num[-1])
    adpt_size = [x for x in range(10, 20000)[::400]]
    print adpt_size, finish_steps, fnl_removed
    fig, ax1 = plt.subplots()
    ax1.bar(adpt_size, fnl_removed, 180,color='y')
    ax1.set_ylabel('removed nodes', color='y')
    ax1.set_xlabel('number of initial adopters')
    ax1.tick_params('y', colors='y')

    ax2 = ax1.twinx()
    ax2.plot(adpt_size, finish_steps, '.-')
    ax2.set_ylabel('step that the infection stops', color='C0')
    ax2.tick_params('y', colors='C0')
    ax2.set_title('SIR - # initial adopters')
    plt.xticks(adpt_size[::6])
    plt.yticks(finish_steps)
    fig.savefig("SIR-adpt_size-finish_step.png")

# largest degree adopters vs. smallest degree adopters
if sys.argv[1] == '5' or sys.argv[1] == 'all':
    max_adopters = []
    min_adopters = []
    degrees = [0 for x in range(number_node)]
    result_in_degree = snap.TIntV()
    result_out_degree = snap.TIntV()
    snap.GetDegSeqV(G, result_in_degree, result_out_degree)
    for i in range(0, result_out_degree.Len()):
        degrees[i] = result_in_degree[i]
    sorted_degrees = sorted(degrees)
    for i in range(init_adopter):
        min_adopters.append(degrees.index(sorted_degrees[i]))
        max_adopters.append(degrees.index(sorted_degrees[-i-1]))

    # ------ adopters with max degrees
    infectious_num, removed_num, final_states = runSIR(G, init_adopter, p, tl, max_adopters)
    steps = [x for x in range (0, len(infectious_num)+1)]
    zero = [0]
    infectious_num = np.array(zero + infectious_num)
    removed_num = np.array(zero + removed_num)
    susceptible_num = np.array([(number_node - infectious_num[x] - removed_num[x]) for x in steps])

    fig, ax = plt.subplots()
    width = 0.2
    p1 = ax.bar(steps, susceptible_num, width)
    p2 = ax.bar(steps, infectious_num, width, bottom=susceptible_num)
    p3 = ax.bar(steps, removed_num, width, bottom=susceptible_num+infectious_num)

    ax.set(xlabel='step number', ylabel='# nodes',
           title='SIR - evolution of the number of nodes in each state')
    plt.xticks(steps)
    ax.legend((p1[0], p2[0], p3[0]), ('susceptible nodes', 'infectious nodes', 'removed nodes'))

    fig.savefig("SIR-step-#state-max_deg.png")

    # ------adopters with min degree
    infectious_num, removed_num, final_states = runSIR(G, init_adopter, p, tl, min_adopters)
    steps = [x for x in range (0, len(infectious_num)+1)]
    zero = [0]
    infectious_num = np.array(zero + infectious_num)
    removed_num = np.array(zero + removed_num)
    susceptible_num = np.array([(number_node - infectious_num[x] - removed_num[x]) for x in steps])

    fig, ax = plt.subplots()
    width = 0.2
    p1 = ax.bar(steps, susceptible_num, width)
    p2 = ax.bar(steps, infectious_num, width, bottom=susceptible_num)
    p3 = ax.bar(steps, removed_num, width, bottom=susceptible_num+infectious_num)

    ax.set(xlabel='step number', ylabel='# nodes',
           title='SIR - evolution of the number of nodes in each state')
    plt.xticks(steps)
    ax.legend((p1[0], p2[0], p3[0]), ('susceptible nodes', 'infectious nodes', 'removed nodes'))

    fig.savefig("SIR-step-#state-min_deg.png")