import sys
import re


graph_file = sys.argv[1]
gaf_file = sys.argv[2]
if len(sys.argv)>3:
	outfile = sys.argv[3]
else:
	outfile = graphfile.split('.gfa')[0]+'-coverage.tsv'

node_lengths, node_coverages = dict(), dict()

with open(graph_file) as f:
	for i,line in enumerate(f):
		if line[0] == 'S':
			name = line.strip().split('\t')[1]
			seq = line.strip().split('\t')[2]			
			node_lengths[name] = len(seq)

with open(gaf_file) as gaf:
	for i,line in enumerate(gaf):
		path = line.strip().split('\t')[5]
		pathlength = int(line.strip().split('\t')[6])
		startposition = int(line.strip().split('\t')[7])
		endposition = int(line.strip().split('\t')[8])

		nodes = re.split('<|>', path)
		nodes = [node for node in nodes if len(node)>= 1]

		#if path consists of one node only:
		if len(nodes) == 1: 
			if nodes[0] not in node_coverages.keys():
				node_coverages[nodes[0]] = float(endposition - startposition)/float(node_lengths[nodes[0]])
			else:
				node_coverages[nodes[0]] += float(endposition - startposition)/float(node_lengths[nodes[0]])
		else:		
			#first node: not all bases are necessarily covered
			firstnode = nodes[0]
			if firstnode not in node_coverages.keys():
				node_coverages[firstnode] = float(node_lengths[firstnode]-startposition)/float(node_lengths[firstnode])
			else:
				node_coverages[firstnode] += float(node_lengths[firstnode]-startposition)/float(node_lengths[firstnode])
			#nodes in the middle of the path: fully covered
			for i in range(1, len(nodes)-1):
				node = nodes[i]
				if node not in node_coverages.keys():
					node_coverages[node] = 1.0
				else:
					node_coverages[node] += 1.0
			#last node: not all bases are necessarily covered
			lastnode = nodes[-1]
			if lastnode not in node_coverages.keys():
				node_coverages[lastnode] = float(node_lengths[lastnode]+endposition-pathlength)/float(node_lengths[lastnode])
			else:
				node_coverages[lastnode] += float(node_lengths[lastnode]+endposition-pathlength)/float(node_lengths[lastnode])

with open(outfile, 'w') as outf:
	for node, length in node_lengths.items():
		if node in node_coverages.keys():
			outf.write(node+'\t'+str(length)+'\t'+str(node_coverages[node])+'\n')
		else:
			outf.write(node+'\t'+str(length)+'\t'+'0.0'+'\n')
