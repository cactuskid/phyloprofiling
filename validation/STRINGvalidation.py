import redis

UNIMAP = True
STRINGMAP = True
data = '/scratch/cluster/monthly/dmoi/'

if UNIMAP:
    # set up redis DB with string data
    r0 = redis.StrictRedis(host='localhost', port=6379, db=0)
    with open(data + 'stringdata/full_uniprot_2_string.04_2015.tsv', 'r')as unimap:
        for i,line in enumerate(unimap):
            try:
                if i > 0:
                    words = line.split('\t')
                    uni = words[1].split('|')[0]
                    string = words[2]
                    r0.set(uni, string)
                if i%10000 == 0:
                    print(i)
            except:
                print(i)

if STRINGMAP:
    r1 = redis.StrictRedis(host='localhost', port=6379, db=1)
    # sort the IDs alphanumerically.
    # protein1 protein2 neighborhood fusion cooccurence coexpression
    # experimental database textmining combined_score
    # 394.NGR_c00010 394.NGR_c33930 0 0 165 0 0 0 145 255
    # 37200000
    # save file line...
    refs = ['neighborhood', 'fusion', 'cooccurence', 'coexpression', 'experimental', 'database', 'textmining', 'combined_score']
    start_line = 0
    with open(data+'/stringdata/protein.links.detailed.v10.5.txt', 'r') as stringAll:
        for i, line in enumerate(stringAll):
            if i > start_line:
                words = line.split()
                IDS = ''.join(sorted([words[0], words[1]]))
                r1.set(IDS, i)
            if i % 1000000 == 0:
                print(i)
