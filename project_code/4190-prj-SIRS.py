import snap
import random
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


def runSIRS (G, init_adopter, p, tl, tr):
    number_node = G.GetNodes()
    number_edge = G.GetEdges()

    initial_list = random.sample(range(number_node), init_adopter)

    # data to return
    infectious_num = []
    removed_num = []
    final_states = []

    # 0: Susceptible, 1:Infectious, 2: Removed
    for NI in G.Nodes():
        nid = NI.GetId()
        state = 0
        if nid in initial_list:
            G.AddIntAttrDatN(nid, 1, "state")
            state = 1
        else:
            G.AddIntAttrDatN(nid, 0, "state")
            state = 0
        G.AddIntAttrDatN(nid, 0, "stepi")
        G.AddIntAttrDatN(nid, 0, "stepr")

    number_R = number_node - len(initial_list)
    number_echo = 0
    ########runs until all nodes are removed########
    ########change this condition to control loop times########
    while number_R != number_node and number_R > 0 and number_echo < 50:
        for NI in G.Nodes():
            nid = NI.GetId()
            ival = G.GetIntAttrDatN(nid, "state")
            # if it's infectious, it will infect its neighbors.
            if ival == 1:
                for nid1 in NI.GetOutEdges():
                    att_value = G.GetIntAttrDatN(nid1, "state")

                    # with probability of p=0.2, a susceptible node will be infected.
                    if (random.randint(0, 4) == 0) and (att_value == 0):
                        # The state 3 is temporary in order to avoid repeated infection.
                        G.AddIntAttrDatN(nid1, 3, "state") # change to 1 later

                # update the step of infection.
                step_this = G.GetIntAttrDatN(nid, "stepi")
                step_this += 1
                G.AddIntAttrDatN(nid, step_this, "stepi")
                if step_this == tl:
                    # The state 4 is temporary in order to avoid repeated infection.
                    G.AddIntAttrDatN(nid, 4, "state") # change to 2 later
                    G.AddIntAttrDatN(nid, 0, "stepi")
            if ival == 2:
                # update the step of removal.
                step_this = G.GetIntAttrDatN(nid, "stepr")
                step_this += 1
                G.AddIntAttrDatN(nid, step_this, "stepr")
                if step_this == tr:
                    # The state 5 is temporary in order to avoid repeated infection.
                    G.AddIntAttrDatN(nid, 5, "state") # change to 0 later
                    G.AddIntAttrDatN(nid, 0, "stepr")

        # Count the number of susceptible or removed nodes.
        number_R_temp = 0
        removed = 0
        for NI in G.Nodes():
            nid = NI.GetId()
            ival = G.GetIntAttrDatN(nid, "state")
            if ival == 3: # The nodes should be updated as infectious now.
                G.AddIntAttrDatN(nid, 1, "state")
                ival = 1
            elif ival == 4:
                G.AddIntAttrDatN(nid, 2, "state")
                ival = 2
            elif ival == 5:
                G.AddIntAttrDatN(nid, 0, "state")
                ival = 0

            if ival == 0 or ival == 2:
                number_R_temp += 1
            if ival == 2:
                removed += 1
        number_R = number_R_temp
        number_echo += 1
        infectious_num.append(number_node - number_R)
        removed_num.append(removed)


    for NI in G.Nodes():
        nid = NI.GetId()
        ival = G.GetIntAttrDatN(nid, "state")
        step_infected = G.GetIntAttrDatN(nid, "stepi")
        step_removed = G.GetIntAttrDatN(nid, "stepr")
        final_states.append((nid, ival, step_infected, step_removed))

    print("One round of simulation completes.")
    return infectious_num, removed_num, final_states



# ---------RUN MODELS & COLLECT DATA----------
if len(sys.argv) != 2:
    print("Usage: <figure #>\n1. tl - infectious node %\n2. tr - infectious node %\n3. p - infectious node %\n4. # nodes in each state - step number\n")

G = snap.LoadEdgeList(snap.PNEANet, "soc-Epinions1.txt", 0, 1, '\t')
number_node = G.GetNodes()

init_adopter = 10
p = 0.2
tl = 2
tr = 2

# tl - infectious %
if sys.argv[1] == '1' or sys.argv[1] == 'all':
    ave_infc_num = [0.05, 0.12, 0.18, 0.23, 0.27]
    # for tl in range(1, 6):
    #     infectious_num, removed_num, final_states = runSIRS(G, init_adopter, p, tl, tr)
    #     last_30_infc = infectious_num[len(infectious_num)-30:len(infectious_num)]
    #     ave_infc_num.append((float(sum(last_30_infc))/30)/number_node)
    infc_interval = [x for x in range(1, 6)]
    print infc_interval, ave_infc_num
    fig, ax = plt.subplots()
    ax.plot(infc_interval, ave_infc_num, '.-')

    ax.set(xlabel='time interval to enter remove state', ylabel='average proportion of infectious nodes',
           title='SIRS - time interval and infectious proportion')
    plt.xticks(infc_interval)
    plt.yticks(ave_infc_num)
    ax.grid()

    fig.savefig("SIRS-tl-infc%.png")

# tr - infectious %
if sys.argv[1] == '2' or sys.argv[1] == 'all':
    ave_infc_num = [0.15,0.12,0.10,0.09,0.08]
    # for tr in range(1, 6):
    #     infectious_num, removed_num, final_states = runSIRS(G, init_adopter, p, tl, tr)
    #     last_30_infc = infectious_num[len(infectious_num)-30:len(infectious_num)]
    #     ave_infc_num.append((float(sum(last_30_infc))/30)/number_node)
    remv_interval = [x for x in range(1, 6)]
    print remv_interval, ave_infc_num
    fig, ax = plt.subplots()
    ax.plot(remv_interval, ave_infc_num, '.-')

    ax.set(xlabel='time interval to leave remove state', ylabel='average proportion of infectious nodes',
           title='SIRS - time interval and infectious proportion')
    plt.xticks(remv_interval)
    plt.yticks(ave_infc_num)
    ax.grid()

    fig.savefig("SIRS-tr-infc%.png")

# p - infectious node %
if sys.argv[1] == '3' or sys.argv[1] == 'all':
    ave_infc_num = [0.11744795441865778, 0.11739919257414216, 0.11768737068227046, 0.11776688323954365, 0.1174519080817266, 0.11820925420735645, 0.11787451073419525, 0.11773261815961376, 0.11770538181402848]
    # for i in range(1, 10):
    #     infectious_num, removed_num, final_states = runSIRS(G, init_adopter, float(i)/10, tl, tr)
    #     last_30_infc = infectious_num[len(infectious_num)-30:len(infectious_num)]
    #     ave_infc_num.append((float(sum(last_30_infc))/30)/number_node)
    contagion_prob  = [float(x)/10 for x in range(1, 10)]
    print contagion_prob, ave_infc_num
    fig, ax = plt.subplots()
    ax.plot(contagion_prob, ave_infc_num, '.-')

    ax.set(xlabel='Contagion Probability', ylabel='average proportion of infectious nodes',
           title='SIRS - contagion probability and infectious proportion')
    plt.xticks(contagion_prob)
    plt.yticks(ave_infc_num)
    ax.grid()

    fig.savefig("SIRS-p-infc%.png")

# infected/removed/susceptible - step number fig
if sys.argv[1] == '4' or sys.argv[1] == 'all':
    infectious_num, removed_num, final_states = runSIRS(G, init_adopter, p, tl, tr)
    steps = [x for x in range (0, len(infectious_num)+1)]
    zero = [0]
    infectious_num = np.array(zero + infectious_num)
    removed_num = np.array(zero + removed_num)
    susceptible_num = np.array([(number_node - infectious_num[x] - removed_num[x]) for x in steps])
    print removed_num

    fig, ax = plt.subplots()
    width = 0.35
    p1 = ax.bar(steps, susceptible_num, width)
    p2 = ax.bar(steps, infectious_num, width, bottom=susceptible_num)
    p3 = ax.bar(steps, removed_num, width, bottom=susceptible_num+infectious_num)

    ax.set(xlabel='step number', ylabel='# nodes',
           title='SIRS - evolution of the number of nodes in each state')
    plt.xticks(steps[0::5])
    ax.legend((p1[0], p2[0], p3[0]), ('susceptible nodes', 'infectious nodes', 'removed nodes'))

    fig.savefig("SIRS-step-#state.png")