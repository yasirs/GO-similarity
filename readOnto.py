import re
import cPickle as pickle
from GOMod import TreeNode

def readOnto(FN="gene_ontology.1_2.obo",verbose=False):
	colon = re.compile(':')
	term = re.compile('\[Term\]')
	GOre = re.compile('GO:[0-9]+')
	FI = open(FN,'r')
	Nodes = []
	Alts = {}
	line = FI.readline()
	while line!='':
		if term.match(line):
			obso=False
			ThisNode = {'parents_id':[],'alts':{}}
			line = FI.readline()
			while (line!='')and(not term.match(line)):
				m = colon.search(line)
				if m:
					tag = line[:m.start()]
					value = line[m.end():].rstrip().lstrip()
					if tag=='id':
						ThisNode['id'] = value
					elif (tag=='is_a')or(tag=='part_of')or(tag=='has_part')or(tag=='regulates')or(tag=='positively_regulates')or(tag=='negatively_regulates'):
						parent = value.split('!')[0].rstrip().lstrip()
						ThisNode.setdefault('parents_id',[]).append(parent)
					elif (tag=='relationship'):
						parent = GOre.search(value).group()
						ThisNode.setdefault('parents_id',[]).append(parent)
					elif (tag=='namespace'):
						ThisNode['namespace'] = value.lstrip()
					elif (tag=='is_obsolete'):
						obso=True
					elif (tag=='alt_id'):
						m = GOre.search(value)
						if m:
							ThisNode['alts'][m.group()]=ThisNode['id']
				else:
					if (line.rstrip()!=''):
						if verbose:
							print "strange line: ",line
				line = FI.readline()
			if ((not obso)and(GOre.match(ThisNode['id']))):
				Nodes.append(ThisNode)
				Alts.update(ThisNode['alts'])
		else:
			if verbose:
				print "bad line : ",line
			line = FI.readline()
	FI.close()
	
	# FO = open('AltGOIds.txt','w')
	# for k in Alts.keys():
	#	FO.write('%s\t%s\n' %(k,Alts[k]))
	# FO.close()
	
	# FO = open('Nodes.txt','w')
	# FO.write(str(len(Nodes))+'\n')
	# for n in Nodes:
	#	FO.write(n['id']+'\n'+n['namespace']+'\n'+str(len(n['parents_id']))+'\n'+'\t'.join(n['parents_id'])+'\n')
	# FO.close()
	
	AllNodes = {}
	Roots = []
	NotLeafs = set()
	for n in Nodes:
		AllNodes[n['id']] = TreeNode(n['id'],n['parents_id'],n['namespace'])
		if len(n['parents_id'])==0:
			Roots.append(n['id'])
		else:
			NotLeafs.update(n['parents_id'])
	Leafs = set(AllNodes.keys()).difference(NotLeafs)
	BadNodes = NotLeafs.difference(set(AllNodes.keys()))
	if verbose:
		print "Number of Leaf Nodes = ",len(Leafs)
		print "Bad Nodes = ", BadNodes
		print "Roots = ", Roots
	
	Onto = {'Nodes':AllNodes,'Roots':Roots}
	
	# FO = open('Onto.pkl','w')
	# pickle.dump(Onto,FO); FO.close()
	return AllNodes, Roots, Alts
