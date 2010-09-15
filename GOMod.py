import cPickle as pickle

class TreeNode:
	def __init__(self,id,parents,namespace):
		self.id=id; self.parents_id=parents; self.namespace=namespace;
		self.Ngenes = 0; self.Ancestors = []; self.doneAnc= False; self.parents_node=[];
	def CompAncestors(self):
		if not self.doneAnc:
			self.Ancestors = list(self.parents_id)
			for p in self.parents_node:
				self.Ancestors.extend(p.getAncestors())
			self.Ancestors.append(self.id)
			self.doneAnc = True
	def getAncestors(self):
		if not self.doneAnc:
			self.CompAncestors()
		return set(self.Ancestors)
	def AddGeneCount(self,g=1):
		self.Ngenes = self.Ngenes + g
		for p in self.parents_node:
			p.AddGeneCount(g)

# def ReadOnto(FN):
#	FI = open(FN,'r')
#	O = pickle.load(FI)
#	FI.close()
#	return O['Nodes'],O['Roots']



		
