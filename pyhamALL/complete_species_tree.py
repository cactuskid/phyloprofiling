import ete3

# load species tree in a ete3 tree object
t = ete3.Tree("speciestree.nwk", format=1, quoted_node_names=True)

# find the node where the node belongs as a child
for node in t.traverse():
    if 'Nanohaloarchaea' in node.name:
        node.add_child(name='Haloredivivus sp. (strain G17)')

# save the corrected tree in a new file
with open(config.working_dir + 'speciestree_corrected.nwk' , 'w') as outfile:
    outfile.write(t.write(format=1,quoted_node_names=True))