import cPickle as pickle
import re
import GOMod
import readOnto
import scipy

class SemanticSimCalculator:

	def __init__(self,ontofile='gene_ontology.1_2.obo',annofile='gene_association.sgd',verbose=False):
		self.ontofile = ontofile
		self.annofile = annofile

		# read the files
		A = readOnto.readOnto(self.ontofile,verbose=verbose)
		#print "type A = ", type(A)
		#print "len A = ", len(A)
		self.Nodes,Roots,AltIDs = A

		if verbose:
			print "Read in the Ontology"
		
		for GOterm in self.Nodes.keys():
			for par in self.Nodes[GOterm].parents_id:
				self.Nodes[GOterm].parents_node.append(self.Nodes[par])



		GOdict = {}
		GOancestor = {}
		for GOterm in self.Nodes.keys():
			GOdict[GOterm] = []
			GOancestor[GOterm] = self.Nodes[GOterm].getAncestors()
		
		if verbose:
			print "Computed Ancestors"
		
		self.GeneAnnot = {}
		if verbose:
			print "reading in annotation and populating the GO tree with annotations"
		# NOTE:: the columns for the GO term and systematic names in the annotation file are defined here
		GOtermind = 4
		systematicNameind = 10
		for line in open(self.annofile,'r'):
		  if (line[0]!='!'):
		     aa=re.split('[\t]',line.rstrip())
		     systnametemp= re.split('\|',aa[systematicNameind])
		     systname=systnametemp[0]
		     GOterm=aa[GOtermind]
		     if (not(GOdict.has_key(GOterm))):
		          if (AltIDs.has_key(GOterm)):
		               GOterm = AltIDs[GOterm]
		     assert(GOdict.has_key(GOterm)), (GOterm + " in the annotation file not recoginzed as a GO term from the ontology!")
		     # if ( not GOdict.has_key(GOterm)):
		     # print "ERROR *******"
		     #	GOdict[GOterm] = []
		     if (systname not in GOdict[GOterm]):
				GOdict[GOterm].append(systname)              # adding systematic name to a GO terms
				self.GeneAnnot.setdefault(systname,{});
				self.GeneAnnot[systname].setdefault(self.Nodes[GOterm].namespace,[]).append(GOterm)
		     for parent in GOancestor[GOterm]:
		             if (not(systname in GOdict[parent])):   # adding systematic to all parent GO terms
		                 GOdict[parent].append(systname)
		
		if verbose:
			print "Read in the Annotation file"
		
		for GOterm in self.Nodes.keys():
			self.Nodes[GOterm].Ngenes = len(GOdict[GOterm])
		if verbose:
			print "done initializing, ready to compute the similarities"


	def getSimGOpair(self,GO1,GO2):
		f = lambda x: (x.namespace) if (type(x)==type(GOMod.TreeNode)) else (self.Nodes[x].namespace)
		assert (f(GO1)==f(GO2)), "GO terms are in different Ontologies! You sucker!!"
		toAnc = lambda x: x.getAncestors() if (type(x)==type(GOMod.TreeNode)) else (self.Nodes[x].getAncestors())
		Anc1 = toAnc(GO1); Anc2 = toAnc(GO2);
		#print "n anc1 = ", len(Anc1), ", n anc2 = ", len(Anc2)
		CommonAnc = list(Anc1.intersection(Anc2))
		t_GO_Ngenes = lambda x: (x.Ngenes,x) if (type(x)==type(GOMod.TreeNode)) else (self.Nodes[x].Ngenes,x)
		Counts = map(t_GO_Ngenes, CommonAnc)
		Counts.sort()
		return (-scipy.log((Counts[0][0]+0.0)/Counts[-1][0]), Counts[0], t_GO_Ngenes(GO1), t_GO_Ngenes(GO2))

	def getSimDictForGenes(self,Gene1,Gene2):
		assert (self.GeneAnnot.has_key(Gene1)), (Gene1 + " not recognized as a gene (need systematic names as given in the annotation file)!")
		assert (self.GeneAnnot.has_key(Gene2)), (Gene2 + " not recognized as a gene (need systematic names as given in the annotation file)!")
		# if (not(self.GeneAnnot.has_key(Gene1) and self.GeneAnnot.has_key(Gene2))):
		#	print "going to return none"
		#	return None
		# print "did not return none"
		ansDict = {'biological_process':(),'molecular_function':(),'cellular_component':()} # n-tuple is (semantic sim,(N genes,LCA), (N genes,GO1), (N genes,GO2))
		for namespace in ansDict.keys():
			pairList = []
			GOList1 = self.GeneAnnot[Gene1][namespace]
			GOList2 = self.GeneAnnot[Gene2][namespace]
			for go1 in GOList1:
				for go2 in GOList2:
					pairList.append(self.getSimGOpair(go1,go2))
			pairList.sort()
			ansDict[namespace] = pairList[-1]
		return ansDict;
	
	def getSimpleSimilarity(self,Gene1,Gene2,namespace='molecular_function'):
	    A=self.getSimDictForGenes(Gene1,Gene2)
	    if A ==None:
	        return None
	    else:
	        return A[namespace][0]


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
		


