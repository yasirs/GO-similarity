import cPickle as pickle
import re
import GOMod


def getMaxSim(GOList1,GOList2,Nodes,namespace='biological_process'):
	f = lambda x: (x.namespace==namespace) if (type(x)==type(GOMod.TreeNode)) else (Nodes[x].namespace==namespace)
	GoodL1 = filter(f ,GOList1)
	GoodL2 = filter(f, GOList2)
	toAnc = lambda x: x.getAncestors() if (type(x)==type(GOMod.TreeNode)) else (Nodes[x].getAncestors())
	Anc1 = set(reduce(lambda x,y: x+y, map(toAnc,GoodL1), []))
	Anc2 = set(reduce(lambda x,y: x+y, map(toAnc,GoodL2), []))
	print "n anc1 = ", len(Anc1), ", n anc2 = ", len(Anc2)
	CommonAnc = Anc1.intersection(Anc2)
	Counts = map(lambda x: x.Ngenes if (type(x)==type(GOMod.TreeNode)) else Nodes[x].Ngenes, CommonAnc)
	print len(Counts), " number of common Ancestors, min = ", min(Counts), ", max = ", max(Counts)
	return (min(Counts)+0.0)/max(Counts)




tairname = re.compile('AT[0-9]*G[0-9]+')

A = pickle.load(open("Onto.pkl",'r'))
Nodes = A['Nodes']
Roots = A['Roots']

for n in Nodes.values():
	ThisParNode = n.parents_node
	for p in n.parents_id:
		ThisParNode.append(Nodes[p])

for n in Nodes.values():
	n.CompAncestors()

Int2Tair = {}
FI = open("arabidopsis.xrefs",'r')
for line in FI:
	if (line[0] != '#'):
		words = line.rstrip().split('\t')
		m = tairname.match(words[7])
		if m:
			Int2Tair[words[1]] = m.group()
FI.close()
Alts = {}
FI=open('AltGOIDs.txt','r')
for line in FI:
	words = line.rstrip().split('\t')
	Alts[words[0]]=words[1]

GeneAnnot = {}
FI = open('gene_association.goa_arabidopsis','r')
for line in FI:
	words = line.rstrip().split('\t')
	Int = words[1]; GOTerm = words[4]
	if words[3]!='NOT':
		if Int2Tair.has_key(Int):
			TairGene = Int2Tair[Int]
			if Nodes.has_key(GOTerm):
				Nodes[GOTerm].AddGeneCount()
			else:
				Nodes[Alts[GOTerm]].AddGeneCount()
			GeneAnnot.setdefault(TairGene,[]).append(GOTerm)

print "similarity = ", getMaxSim(GeneAnnot["At1g59750".upper()],GeneAnnot["At2g46530".upper()],Nodes)


		





