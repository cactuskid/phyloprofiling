def Tree2Hashes(treemap, fam=None, LSH=None):

    eventdict = tree2eventdict(treemap)



        if LSH is not None:
            LSH.insert(lminHashName, lminHash)


        buf = bytearray(lminHash.bytesize())
        lminHash.serialize(buf)

        hashes[array] = buf


    for j in range(1,len(eventdict.keys())):
        for i in itertools.combinations(eventdict.keys(), j+1):
            combName = str(fam)
            minHash = datasketch.MinHash(num_perm=128)
            for array in i:
                combName += '-'+array
                minHash.merge(hashesDict[array])

            lminHashDict[combName] = minHash
            lminHash = datasketch.LeanMinHash(minHash)
        if LSH:
            LSH.insert( combName , lminHash)

    return hashes , lminHashDict