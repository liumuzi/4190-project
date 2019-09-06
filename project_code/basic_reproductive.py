import snap
import random
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


def runSIR (G, init_adopter, p, tl):
    number_node = G.GetNodes()
    number_edge = G.GetEdges()

    initial_list = random.sample(range(number_node), init_adopter)

    # data to return
    infectious_num = []
    removed_num = []
    final_states = []
    numberlist = []
    new_I_list = []

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
    while number_echo < 50:
        number_new_I = 0
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
                        number_new_I += 1
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
        number_I = 0
        for NI in G.Nodes():
            nid = NI.GetId()
            ival = G.GetIntAttrDatN(nid, "state")
            infected = G.GetIntAttrDatN(nid, "infected")
            if ival == 3: # The nodes should be updated as infectious now.
                G.AddIntAttrDatN(nid, 1, "state")
                ival = 1
            elif ival == 4:
                G.AddIntAttrDatN(nid, 0, "state")
                ival = 0          
            if ival == 0:
                number_R_temp += 1
            if infected == 1:
                number_I += 1
        number_R = number_R_temp
        number_echo += 1
        infectious_num.append(number_node - number_R)
        numberlist.append(number_I)
        new_I_list.append(number_new_I)
    print("One round of simulation completes.")
    return numberlist,new_I_list

'''
    for NI in G.Nodes():
        nid = NI.GetId()
        ival = G.GetIntAttrDatN(nid, "state")
        step_number = G.GetIntAttrDatN(nid, "step")
        final_states.append((nid, ival, step_number))
        # print(nid, ival, step_number)
'''


# ---------RUN MODELS & COLLECT DATA----------
if len(sys.argv) != 2:
    print("Add 1")

G = snap.LoadEdgeList(snap.PNEANet, "soc-Epinions1.txt", 0, 1, '\t')
number_node = G.GetNodes()

init_adopter = 1
p = 0.15
tl = 1

if sys.argv[1] == '1':
    steps = [x for x in range(50)]
    #t_total_list = [0 for x in range(50)]
    t_new_list = [0 for x in range(50)]
    for i in range(10):
        #temp_total_list = t_total_list
        temp_new_list = t_new_list
        total_infect_list, new_infect_list = runSIR(G,init_adopter,p,tl)
        print(new_infect_list)
        #t_total_list = [temp_total_list[k]+total_infect_list[k] for k in range(50)]
        t_new_list = [temp_new_list[k]+new_infect_list[k] for k in range(50)]
    #f_total_list = [float(t_total_list[k])/10.0 for k in range(50)]
    f_new_list = [float(t_new_list[k])/10.0 for k in range(50)]
    print(f_new_list)
    fig, ax = plt.subplots()
    ax.plot(steps,f_new_list,color='royalblue',marker='o')
    ax.set(xlabel='step number', ylabel='# nodes',title='Newly infected nodes vs. Steps (p = 0.15)')
    plt.xticks([0,5,10,15,20,25,30,35,40,45,50])
    fig.savefig("basic_reproductive_0.15.png")
