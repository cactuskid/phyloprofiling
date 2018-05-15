def DFTree2Hashes(row):

    fam, treemap = row.tolist()

    eventdict = { 'presence':[] , 'gain':[] , 'loss':[] , 'duplication':[]}
    if treemap is not None:
        for node in treemap.traverse():
        # traverse() returns an iterator to traverse the tree structure
        # strategy:"levelorder" by default; nodes are visited in order from root to leaves
        # it return treeNode instances
            if not node.is_root():
                if node.nbr_genes >0:
                    eventdict['presence'].append('P'+node.name)
                if node.dupl > 0:
                    eventdict['duplication'].append('D'+node.name)
                if node.lost > 0:
                    eventdict['loss'].append('L'+node.name)
            else:
                eventdict['gain'].append('G'+node.name)

        hashes = {}

        hashesDict = {}

        lminHashDict = {}

        for array in eventdict:
            eventdict[array] = set(eventdict[array])

            minHash = datasketch.MinHash(num_perm=128)

            for element in eventdict[array]:

                minHash.update(element.encode())


            hashesDict[array] = minHash

            lminHash = datasketch.LeanMinHash(minHash)

            lminHashName = str(fam)+'-'+array

            lminHashDict[lminHashName] = lminHash

            #lminHash.bytesize() --> computes the bitesize
            #buf = bytearray(lminHash.bytesize())
            #lminHash.serialize(buf)

            #hashes[array] = buf

            hashes[array] = lminHash.hashvalues


			buf = bytearray(lminHash.bytesize())
			lminHash.serialize(buf)

			hashes[array] = buf
        #hashmat = np.hstack(buffers)

        for j in range(1,len(eventdict.keys())):
            for i in itertools.combinations(eventdict.keys(), j+1):
                combName = str(fam)
                minHash = datasketch.MinHash(num_perm=128)
                for array in i:
                    combName += '-'+array
                    minHash.merge(hashesDict[array])

                lminHashDict[combName] = minHash
                lminHash = datasketch.LeanMinHash(minHash)

        return {'hashes':hashes , 'dict': lminHashDict}
    else:
        return None