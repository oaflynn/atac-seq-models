import matplotlib.pyplot as plt
import pandas as pd
import math
import numpy as np

# Input
NARROWPEAK_NAME = 'narrowPeak/72hr_rep1.narrowPeak'
OUT_NAME = 'Openness output/72hr_rep1_peaks'
L_0 = 132000000 # Genome size
PSEUDOCOUNT = 5 # delta term

# Go through each entry in narrowPeak file
Y = 0 # Initialize genome peak aggregate
peak_data = [] # List to store (peak height, peak length, p-value, chrom, start) tuples
narrowPeak_file = open(NARROWPEAK_NAME, 'r')
for peak in narrowPeak_file:
    data = peak.split()
    height = float(data[6])
    length = float(data[2]) - float(data[1])
    Y += height # Construct Y term

    pvalue = float(data[7])
    chrom = data[0]
    start = float(data[2])
    peak_data.append((height, length, pvalue, chrom, start))
narrowPeak_file.close()

# Compute openness score for each peak
o_scores = []
csv = open(OUT_NAME+'.csv', 'w')
csv.write('chrom,start,end,height,length,pvalue,openness\n')
for height, length, pvalue, chrom, start in peak_data:
    X = height * (1-(10**(-1*pvalue))) # Incorproate peak p-value
    L = length
    o_score = (X / L) / ((Y+PSEUDOCOUNT) / L_0)
    o_scores.append(o_score)

    csv.write(chrom+','+str(start)+','+str(start+length-1)+','+str(height)+','+str(length)+','+str(pvalue)+','+str(o_score)+'\n')
csv.close()

# Some descriptive statistics
print('Descriptive statistics for ' + OUT_NAME)
print('Min openness score: ', min(o_scores))
print('Max openness score: ', max(o_scores))
print('Avg openness score: ', sum(o_scores)/len(o_scores))
print('Median openness score: ', np.median(o_scores))
print('Standard deviation: ', np.std(o_scores))

# Plot histogram
plt.hist(o_scores, bins='auto')
plt.xlabel('Openness score (log10)')
plt.ylabel('# of peaks')
plt.figtext(0.99, 0.01, str(len(o_scores)) + ' peaks', horizontalalignment='right')
plt.savefig(OUT_NAME+'.png')