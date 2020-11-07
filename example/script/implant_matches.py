#r_map = '/data6/playyard/ukbiobank/tmp/debug_rapid_misc/decode_maps_hg19_filtered/decode_22_hg19.txt'
r_map = 'decode_22_hg19.txt'
num_matches = 100
min_length = 5
min_breaking = 2  # currently only the middle point is the breaking point

'chr22   17582788        17582789        0.013987185309072841    0.21706195005459356'
f = open(r_map)
sites = []
site_pos = {}
for line in f:
    vals = line.split()
    pos = int(vals[1])
    pos_g = float(vals[4])
    sites.append(pos)
    site_pos[pos] = pos_g
f.close()

import random
panel = []
num_samples = 200
for n in range(0, num_samples):
    vals = ['1|1'] * len(sites)
    for s in range(0, len(sites)):
        vals[s] = random.choice(['0|0', '0|1', '1|1', '1|0'])
    panel.append(vals)

start_end_matches = []

while (len(start_end_matches) < num_matches):
    sp = random.randint(0, len(site_pos))
    gp = site_pos[sites[sp]]
    dist = 0
    ge = sp
    while (dist < min_length):
        dist = site_pos[sites[ge]] - site_pos[sites[sp]]
        ge = ge + 1
        if (ge > len(sites) - 1): break
    if (dist >= min_length):
        start_end_matches.append([sp, ge])

individuals = set()
overlap_cases = 4
gap_cases = 5
overlap_count = 0
#overlap = 3
gap = 3
for se in start_end_matches:
    s1 = random.randint(0, num_samples - 1)
    s2 = random.randint(0, num_samples - 1)

    while (s1 == s2 or s1 in individuals or s2 in individuals):
        s2 = random.randint(0, num_samples - 1)
        s1 = random.randint(0, num_samples - 1)

    individuals.add(s1)
    individuals.add(s2)
    midpoint = int((se[0] + se[1]) / 2)
    vals = ['0'] * len(sites)
    for s in range(se[0], midpoint):
        sr = random.choice(['0', '1'])
        r1 = random.choice(['0', '1'])
        r2 = random.choice(['0', '1'])

        panel[s1][s] = sr + "|" + r1
        panel[s2][s] = sr + "|" + r2

    # inserting strict boundaries
    # make sure the match ends
    panel[s1][midpoint] = '0' + panel[s1][midpoint][1:]
    panel[s2][midpoint] = '1' + panel[s1][midpoint][1:]

    # NoGap+overlap Case
    # make sure the match starts exactly from a given location
    # allowing some overlap
    if overlap_count < 4:
        #for s in range(midpoint - overlap, se[1]):
        for s in range(midpoint + gap, se[1] + 5):
            sr = random.choice(['0', '1'])

            panel[s1][s] = panel[s1][s][:2] + sr
            panel[s2][s] = panel[s2][s][:2] + sr

        # insert mismatch
        #panel[s1][midpoint - overlap -
        #          1] = panel[s1][midpoint - overlap - 1][:2] + '0'
        #panel[s2][midpoint - overlap -
        #          1] = panel[s2][midpoint - overlap - 1][:2] + '1'

        panel[s1][midpoint + gap - 1] = panel[s1][midpoint + gap - 1][:2] + '0'
        panel[s2][midpoint + gap - 1] = panel[s2][midpoint + gap - 1][:2] + '1'
        print('NA' + str(s1), 'NA' + str(s2), se[0], midpoint + gap, se[1] + 5,
              sites[se[0]], sites[midpoint + gap], sites[se[1] + 5])
        overlap_count += 1
    else:
        for s in range(midpoint, se[1]):
            sr = random.choice(['0', '1'])
            r1 = random.choice(['0', '1'])
            r2 = random.choice(['0', '1'])

            panel[s1][s] = r1 + "|" + sr
            panel[s2][s] = r2 + "|" + sr

        # make sure the match ends
        panel[s1][se[1]] = r1 + "|" + '0'
        panel[s2][se[1]] = r2 + "|" + '1'

        # NoGap Case
        # make sure the match starts exactly from a given location
        panel[s1][midpoint - 1] = panel[s1][midpoint - 1][:2] + '0'
        panel[s2][midpoint - 1] = panel[s2][midpoint - 1][:2] + '1'

        print('NA' + str(s1), 'NA' + str(s2), se[0], midpoint, se[1],
              sites[se[0]], sites[midpoint], sites[se[1]])
    #print (panel[s1][se[0]:se[1]])
    #print (panel[s2][se[0]:se[1]])

#output_file = 'test_chr22woverlap.vcf'
output_file = 'test_chr22wgap.vcf'
f_output = open(output_file, 'w+')
f_output.write('##fileformat=VCFv4.1\n')
f_output.write('##source=pseq\n')
f_output.write('##FILTER=<ID=PASS,Description="Passed variant FILTERs">\n')
f_output.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT')

for i in range(0, num_samples):
    f_output.write('\t' + 'NA' + str(i))
f_output.write('\n')

counter = 0
for s in range(0, len(sites)):
    f_output.write('chr1' + '\t' + str(sites[s]) + '\t' + 'rs' + str(counter) +
                   '\t' + 'A' + '\t' + 'G' + '\t'
                   '.' + '\t' + 'PASS' + '\t' + '.' + '\t' + 'GT\t')
    for i in range(0, num_samples):
        f_output.write(panel[i][s] + "\t")
    f_output.write("\n")
    counter = counter + 1

f_output.close()
