import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.stats as stats
import matplotlib.image as mpimg
import sys
import seaborn as sns

import matplotlib.patches as mplpatches

try:
    psl = open(sys.argv[1])
    numreads = 50
except:
    sys.stderr.write('python script.py psl out.png\n')
    sys.exit(1)

try:
    plt.style.use('~/bme/BME163')
except:
    sys.stderr.write('cannot find stylesheet\n')
    pass

def get_junctions(line):
    junctions = set()
    starts = [int(n) for n in line[20].split(',')[:-1]]
    sizes = [int(n) for n in line[18].split(',')[:-1]]
    if len(starts) == 1:
        return
    for b in range(len(starts)-1):
        junctions.add((starts[b]+sizes[b], starts[b+1]))
    return junctions

def parse_psl(psl, names=False, dontadd=set()):
    info = []
    lowbound, upbound = 1e9, 0
    chrom = ''
    isonames = []
    for line in psl:
        if len(info) > numreads:  # max number of reads to plot
            break
        line = line.rstrip().split('\t')
        if line[9] in dontadd:
            continue
        lowbound = min(lowbound, int(line[15]))
        upbound = max(upbound, int(line[16]))
        blocksizes = [int(n) for n in line[18].split(',')[:-1]]
        blockstarts = [int(n) for n in line[20].split(',')[:-1]]
        info += [[blocksizes, blockstarts]]
        strand = line[8]
        isonames += [line[9]]
    upbound += 100
    lowbound -= 100
    for i in range(len(info)):
        info[i][1] = [n - lowbound for n in info[i][1]]
    if names:
        return info, isonames
    else:
        return info, lowbound, upbound, strand, isonames

def parse_junctions(psl, names=False, plotany=True):
    info = []
    usednames = []
    for line in psl:
        if len(info) > numreads:
            break
        line = line.rstrip().split('\t')
        if line[13] != chrom or int(line[15]) > upper or int(line[16]) < lower:
            continue
        blocksizes = [int(n) for n in line[18].split(',')[:-1]]
        blockstarts = [int(n) - lower for n in line[20].split(',')[:-1]]

        prevEnd = blockstarts[0] + blocksizes[0]

        for size, start in zip(blocksizes[1:], blockstarts[1:]):
            junc = [[100]*2, [prevEnd-100, start]]
            if junc not in info:
                info += [junc]
            prevEnd = start + size
    return info

def parse_gtf(gtf):
    info = []
    chrom, lower, upper ='chr1', 206904300, 206922100  # fcmr
    chrom, lower, upper = 'chr15', 44711477, 44718170  # b2m
    chrom, lower, upper = 'chr11', 61964550, 61967660  # fth1
    current_transcript = ''
    for line in gtf:
        line = line.rstrip().split('\t')
        if len(line) < 8:
            continue
        start = line[8].find('transcript_id')
        if start < 0:  # not a transcript
            continue
        transcript_id = line[8][start+len('transcript_id')+2:]
        transcript_id = transcript_id[:transcript_id.find('";')]
        if transcript_id != current_transcript:
            current_transcript = transcript_id
            info += [[[],[],[]]]  # size_list, start_list, utr_boolean
        if line[2] == 'exon':
            info[-1][2] += [0]  # boolean
        elif line[2] == 'UTR':
            info[-1][2] += [1]  # UTRs are plotted with smaller blockheights
        else:
            continue
        info[-1][0] += [int(line[4]) - int(line[3])]
        info[-1][1] += [int(line[3]) - lower]
    return info, lower, upper, chrom

def pack(data, rev=True, color=False, tosort = True):
    starts = [max(d[1]) for d in data] if rev else [min(d[1]) for d in data] # sort by right or left end
    if tosort:
        data = [d for (s,d) in sorted(zip(starts, data))]
    else:
        data = [d for s, d in zip(starts, data)]
    packed = [[ data[0] ]]
    ends = [max(data[0][1]) + data[0][0][data[0][1].index(max(data[0][1]))]]
    for i in range(1, len(data)):
        min_start = min(data[i][1])
        end = max(data[i][1]) + data[i][0][data[i][1].index(max(data[i][1]))]
        pos = -1
        for j in range(len(packed)):
            if ends[j] + 5 < min_start:  # added the plus 5 for spacing between reads
                pos = j  # pack this read with the read at position j
                break
        if pos >= 0:
            packed[pos] += [data[i]]
            ends[pos] = end
        else:
            packed += [[data[i]]]
            ends += [end]
    return packed

def plot_blocks(data, panel, names, utr=False, height=.5, l=0.8, color=False):
    panel.set_xlim(1, upper - lower)
    panel.set_ylim(-1, len(data) * 2)
    # panel.set_ylim(-1, 10)

    panel.tick_params(axis='both', which='both',\
                       bottom=False, labelbottom=False,\
                       left=False, labelleft=False,\
                       right=False, labelright=False,\
                       top=False, labeltop=False)

    di = 0  # data index
    ni = 0  # nameindex
    for i in range(len(data)*2):  # each line
        if i % 2 == 0:
            panel.text(data[di][0][1][0], i-height/2 + 1, names[ni], fontsize=6, ha='left', va='center')
            ni += 1
            continue
        read = data[di]
        di += 1
        for j in range(len(read)):  # each read on each line
            line = read[j]
            sizes = line[0]
            starts = line[1]
            panel.plot([min(starts)+2, max(starts)], [i - 1]*2, 'k-', lw=l)
            for k in range(len(sizes)):  # each block of each read
                if utr and line[2][k]:
                    rectangle1 = mplpatches.Rectangle([starts[k], i-height/2], \
                        sizes[k], height, facecolor='white', linewidth=0, zorder=11)
                    rectangle2 = mplpatches.Rectangle([starts[k], i-.25/2],\
                        sizes[k], 0.25, facecolor='black', linewidth=0, zorder=12)
                    panel.add_patch(rectangle1)
                    panel.add_patch(rectangle2)
                elif color:
                    rectangle = mplpatches.Rectangle([starts[k], i - height/2 - 1], \
                        sizes[k], height, facecolor=color[di], linewidth=0, zorder=10)
                    panel.add_patch(rectangle)
                    print(color[di])             
                else:
                    rectangle = mplpatches.Rectangle([starts[k], i - height/2], \
                        sizes[k], height, facecolor='black', linewidth=0)
                    panel.add_patch(rectangle)

purple, lightblue, orange, red = "#9b59b6", "#3498db", "#ffab40", "#d3494e"
grayblue, darkblue, yellow = "#95a5a6", "#34495e", "#ffd42a"
flatui_mod = sns.color_palette([darkblue,red,yellow])
colors = ["windows blue", "faded green", 'amber', 'dusty purple']
colors = sns.xkcd_palette(colors) + \
        flatui_mod


# bars
fig_0 = plt.figure(figsize=(10, 2))
panel0 = plt.axes([0.005, 0.015, .99, 0.97], frameon=True)  # annotation

reads, lower, upper, strand, names = parse_psl(psl)
upper += 1
if strand == '-':
    print('- strand, inverted')
    panel0.invert_xaxis()
elif strand != '+':
    print('ambiguous strand assumed to be +')

plot_blocks(pack(reads,tosort=False), panel0, names, l=0.4, color=colors)  # plot as is


plt.savefig(sys.argv[1][:-4]+'_bars.pdf', transparent=True)
