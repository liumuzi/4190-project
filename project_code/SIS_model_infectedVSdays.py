import snap
import random
import numpy as np
import matplotlib.pyplot as plt

'''
number_node = 10
number_edge = 30
G = snap.GenRndGnm(snap.PNEANet, number_node, number_edge)
'''
G = snap.LoadEdgeList(snap.PNEANet,"./soc-Epinions1.txt",0,1)
number_node = G.GetNodes()
number_edge = G.GetEdges()

#initial_list = [1, 2, 3]
initial_list = [random.randint(0,number_node-1)] # random select one as initial adopter
#initial_list = random.sample(range(number_node), 10)

p = 0.2 # 0.2 0.4 0.6 0.8
daylist = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
numberlist = [1]
tl = 2

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

    string_node = str(NI.GetId()) + ';' + str(0) + ';' + str(state)
    labels[NI.GetId()] = string_node


#snap.DrawGViz(G, snap.gvlDot, "temp_graph.png", "graph 1", labels)
number_R = number_node - len(initial_list)
number_echo = 1 # number of rounds (complete days?)
#while number_R != number_node:
while number_echo < 17:
    for NI in G.Nodes():
        nid = NI.GetId()
        ival = G.GetIntAttrDatN(nid, "state")
        # if it's infectious, it will infect its neighbors.
        if ival == 1:
            for nid1 in NI.GetOutEdges():
                att_value = G.GetIntAttrDatN(nid1, "state")

                # with probability of p, a susceptible node will be infected.
                if (np.random.choice(np.arange(0,2), p=[p, 1-p]) == 0) and (att_value == 0):
                    # The state 3 is temporary in order to avoid repeated infection.
                    G.AddIntAttrDatN(nid1, 3, "state")
                    G.AddIntAttrDatN(nid1, 1, "infected")

            # update the step of infection.
            step_this = G.GetIntAttrDatN(nid, "step")
            step_this += 1
            if step_this == tl:
                # The state 4 is temporary in order to avoid repeated infection.
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

        # update the label for each node for visualization.
        step_number = G.GetIntAttrDatN(nid, "step")
        string_node = str(NI.GetId()) + ';' + str(step_number) + ';' + str(ival)
        labels[NI.GetId()] = string_node
        #snap.DrawGViz(G, snap.gvlDot, "temp_graph%d.png" % number_echo, "graph %d" % number_echo, labels)

        if ival == 0:
            number_R_temp += 1
    number_R = number_R_temp
    number_echo += 1
    #print(number_R)

    # Count the number of people that have been infected
    number_I = 0
    for NI in G.Nodes():
        nid = NI.GetId()
        infected = G.GetIntAttrDatN(nid, "infected")
        if infected == 1:
            number_I += 1
    numberlist.append(number_I)
print(numberlist)

plt.plot(daylist,numberlist,marker='o')
plt.xlabel("Lasting Days")
plt.ylabel("Number of Infected People")
plt.title("Number of Infected People vs Lasting Days - SIS model (p=0.2)")
plt.savefig("infected_days_SIS_0.2.jpg")