import snap
import random

G = snap.LoadEdgeList(snap.PNEANet,"./soc-Epinions1.txt",0,1)
number_node = G.GetNodes()
number_edge = G.GetEdges()

total_outdeg = 0
for NI in G.Nodes():
	total_outdeg += NI.GetOutDeg()
avg_outdeg = float(total_outdeg) / float(number_node)