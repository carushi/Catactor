with open('scobj_list.txt') as f:
    for line in sorted(f.readlines()):
        # print(line)
        contents = line.split('_')
        tdidf = ('tfidf' if 'tdidf' in line else '')
        if 'bin_bin' in line:
            continue
        binary = ('binary' if 'bin_scanpy' in line else '')
        print(contents[0], contents[1], tdidf, binary, line.rstrip('\n'))