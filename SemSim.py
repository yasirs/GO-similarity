import cPickle as pickle
import re
import GOMod
import readOnto
import scipy




A = readOnto.readOnto('gene_ontology.1_2.obo')
print "type A = ", type(A)
print "len A = ", len(A)
Nodes,Roots,AltIDs = A

print "Read in the Ontology"


for GOterm in Nodes.keys():
	for par in Nodes[GOterm].parents_id:
		Nodes[GOterm].parents_node.append(Nodes[par])


GOdict = {}
GOancestor = {}
for GOterm in Nodes.keys():
	GOdict[GOterm] = []
	GOancestor[GOterm] = Nodes[GOterm].getAncestors()

print "Computed Ancestors"

GeneAnnot = {}


count=1;GOtermind=4;
for line in open('gene_association.sgd','r'):
  if (count>27):
     aa=re.split('[\t]',line.rstrip())
     ORFnametemp= re.split('\|',aa[10])
     ORFname=ORFnametemp[0]
     GOterm=aa[GOtermind]
     if (not(GOdict.has_key(GOterm))):
          if (AltIDs.has_key(GOterm)):
               GOterm = AltIDs[GOterm]
     if ( not GOdict.has_key(GOterm)):
	print "ERROR *******"
	GOdict[GOterm] = []
     if (ORFname not in GOdict[GOterm]):
		GOdict[GOterm].append(ORFname)              # adding ORF to a GO terms
		GeneAnnot.setdefault(ORFname,{});
		GeneAnnot[ORFname].setdefault(Nodes[GOterm].namespace,[]).append(GOterm)
     for parent in GOancestor[GOterm]:
             if (not(ORFname in GOdict[parent])):   # adding ORF to all parent GO terms
                 GOdict[parent].append(ORFname)
  count=count+1

print "Read in the Annotation file"

for GOterm in Nodes.keys():
	Nodes[GOterm].Ngenes = len(GOdict[GOterm])


def getSimGOpair(GO1,GO2):
	f = lambda x: (x.namespace) if (type(x)==type(GOMod.TreeNode)) else (Nodes[x].namespace)
	assert (f(GO1)==f(GO2)), "GO terms are in different Ontologies! You sucker!!"
	toAnc = lambda x: x.getAncestors() if (type(x)==type(GOMod.TreeNode)) else (Nodes[x].getAncestors())
	Anc1 = toAnc(GO1); Anc2 = toAnc(GO2);
	#print "n anc1 = ", len(Anc1), ", n anc2 = ", len(Anc2)
	CommonAnc = list(Anc1.intersection(Anc2))
	t_GO_Ngenes = lambda x: (x.Ngenes,x) if (type(x)==type(GOMod.TreeNode)) else (Nodes[x].Ngenes,x)
	Counts = map(t_GO_Ngenes, CommonAnc)
	Counts.sort()
	return (-scipy.log((Counts[0][0]+0.0)/Counts[-1][0]), Counts[0], t_GO_Ngenes(GO1), t_GO_Ngenes(GO2))



#def getMaxSim(GOList1,GOList2,namespace='biological_process'):
#	f = lambda x: (x.namespace==namespace) if (type(x)==type(GOMod.TreeNode)) else (Nodes[x].namespace==namespace)
#	GoodL1 = filter(f ,GOList1)
#	GoodL2 = filter(f, GOList2)
#	toAnc = lambda x: x.getAncestors() if (type(x)==type(GOMod.TreeNode)) else (Nodes[x].getAncestors())
#	Anc1 = set(reduce(lambda x,y: x+y, map(toAnc,GoodL1), []))
#	Anc2 = set(reduce(lambda x,y: x+y, map(toAnc,GoodL2), []))
#	print "n anc1 = ", len(Anc1), ", n anc2 = ", len(Anc2)
#	CommonAnc = Anc1.intersection(Anc2)
#	Counts = map(lambda x: x.Ngenes if (type(x)==type(GOMod.TreeNode)) else Nodes[x].Ngenes, CommonAnc)
#	print len(Counts), " number of common Ancestors, min = ", min(Counts), ", max = ", max(Counts)
#	return (min(Counts)+0.0)/max(Counts)



#def getSemSim(GO1,GO2):
#	getMaxSim([GO1],[GO2],Nodes)
		
def getSimDictForGenes(Gene1,Gene2):
	if (not(GeneAnnot.has_key(Gene1) and GeneAnnot.has_key(Gene2))):
		print "going to return none"
		return None
	print "did not return none"
	ansDict = {'biological_process':(),'molecular_function':(),'cellular_component':()} # n-tuple is (semantic sim,(N genes,LCA), (N genes,GO1), (N genes,GO2))
	for namespace in ansDict.keys():
		pairList = []
		GOList1 = GeneAnnot[Gene1][namespace]
		GOList2 = GeneAnnot[Gene2][namespace]
		for go1 in GOList1:
			for go2 in GOList2:
				pairList.append(getSimGOpair(go1,go2))
		pairList.sort()
		ansDict[namespace] = pairList[-1]
	return ansDict;

def getSimpleSimilarity(Gene1,Gene2,namespace='molecular_function'):
    A=getSimDictForGenes(Gene1,Gene2)
    if A ==None:
        return None
    else:
        return A[namespace][0]


