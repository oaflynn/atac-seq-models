import matplotlib.pyplot as plt
import pandas as pd
import math
from bioinfokit import analys, visuz

# Input
NARROWPEAK_NAME = 'narrowPeak/48hr_rep1.narrowPeak'
GTF_NAME = 'dm6.ensGene.gtf'
OUT_NAME = 'Openness output/500promoter500gene_48hr_rep1' # Output path WITHOUT extension
L_0 = 132000000 # Genome size
PSEUDOCOUNT = 5 # delta term
PROMOTER_BUFFER = 500 # The openness score will be taken from the region PROMOTER_BUFFER base pairs before and after the start of the gene

PLOT_OPENNESS_DISTS = True
PLOT_PVAL_DIST = True
PLOT_ACCESSIBILITY_DISTS = True
VOLCANO_PLOT = True # Fold change of openness from including vs excluding promoter

# A dictionary structure to keep track of all of the peaks by their location
# {chrom: {start_bp: entry_data}}
chrom_dict = {}

# Go through each entry in narrowPeak file
Y = 0 # Initialize genome peak aggregate
narrowPeak_file = open(NARROWPEAK_NAME, 'r')
for peak in narrowPeak_file:
    data = peak.split()
    Y += float(data[6]) # Construct Y term
    try: # May need to initialize an empty dictionary if an entry for this chromosome doesn't exist yet
        chrom_dict[data[0]][data[1]] = data # Add the peak to the dictionary
    except KeyError:
        chrom_dict[data[0]] = {}
        chrom_dict[data[0]][data[1]] = data
narrowPeak_file.close()

# Prepare output file
out_file = open(OUT_NAME+'.csv', 'w')
out_file.write('chrom,gene_id,openness_score,openness_gene,openness_promoter,avg_pValue\n')

# Iterate through gtf file to match peaks to genes
gtf_file = open(GTF_NAME, 'r')
openness_hist_data = []
openness_gene_hist_data = []
openness_promoter_hist_data = []
p_hist_data = []
volcano_df = pd.DataFrame(columns=['GeneNames', 'value1', 'value2', 'log2fc', 'p-value'])
no_peaks = 0 # Keep track of how many genes have no peaks
for entry in gtf_file:
    # Only look at entries describing genes, which are marked 'transcript'
    if 'transcript\t' not in entry:
        continue

    # Prepare gene data
    data = entry.split()
    gene_start = int(data[3])
    gene_end = int(data[4])

    # Find peaks relavent to this gene
    try: # If the chromosome had no narrowPeak data, record the openness score as 0
        peaks_dict = chrom_dict[data[0]] # Peaks on the chromosome
    except KeyError:
        line = data[0] + ',' + data[9][1:-2] + ',0\n'
        out_file.write(line)
        continue
    front_buffer = PROMOTER_BUFFER if data[6] == '+' else 0 # The strand direction shows which side of the gene the promoter will be found on
    back_buffer = PROMOTER_BUFFER if data[6] == '-' else 0
    gene_peaks = [peaks_dict[key] for key in peaks_dict.keys() if (int(key) >= gene_start and int(key) < gene_start+PROMOTER_BUFFER and front_buffer > 0) or (int(key) >= gene_end-PROMOTER_BUFFER and int(key) < gene_end and back_buffer > 0)] # Peaks starting within gene body region 
    promoter_peaks = [peaks_dict[key] for key in peaks_dict.keys() if (int(key) >= gene_start-front_buffer and int(key) < gene_start) or (int(key) < gene_end+back_buffer and int(key) >= gene_end)] # Peaks starting within the promoter region
    all_peaks = gene_peaks + promoter_peaks

    # Process peaks within gene area
    X = sum([float(peak_data[6])*(1-10**(-1*float(peak_data[7]))) for peak_data in all_peaks]) # Aggregate gene peaks, factoring in p-value
    X_gene = sum([float(peak_data[6])*(1-10**(-1*float(peak_data[7]))) for peak_data in gene_peaks])
    X_promoter = sum([float(peak_data[6])*(1-10**(-1*float(peak_data[7]))) for peak_data in promoter_peaks])
    avg_p = [float(peak_data[7]) for peak_data in all_peaks] # Average p-value for peaks
    if len(avg_p) > 0:
        avg_p = sum(avg_p) / len(avg_p)
    else:
        avg_p = 0

    # Compute openness
    L = gene_end - gene_start # Gene length
    openness_score = (X / (L+PROMOTER_BUFFER)) / ((Y + PSEUDOCOUNT) / L_0)
    openness_gene = (X_gene / L) / ((Y + PSEUDOCOUNT) / L_0)
    openness_promoter = (X_promoter / PROMOTER_BUFFER) / ((Y + PSEUDOCOUNT) / L_0)

    # Write to file
    line = data[0] + ',' + data[9][1:-2] + ',' + str(openness_score) + ',' + str(openness_gene) + ',' + str(openness_promoter) + ',' + str(avg_p) + '\n'
    out_file.write(line)

    # Save openness for histogram
    if openness_score > 0:
        openness_hist_data.append(math.log10(openness_score))
        p_hist_data.append(avg_p)
    else:
        no_peaks += 1
    if openness_gene > 0:
        openness_gene_hist_data.append(math.log10(openness_gene))
    if openness_promoter > 0:
        openness_promoter_hist_data.append(math.log10(openness_promoter))

    # Save data for volcano plot, fold change from only gene body to gene+promoter
    if openness_gene > 0:
        fold_change = openness_score / openness_gene
        row = {'GeneNames':data[9][1:-2], 'value1':openness_score, 'value2':openness_gene, 'log2fc':math.log2(fold_change), 'p-value':10**(-1*avg_p)}
        volcano_df = volcano_df.append(row, ignore_index=True)

out_file.close()
gtf_file.close()

# Append openness probabilities
df = pd.read_csv(OUT_NAME+'.csv')
max_openness = max(df['openness_score'])
sum_openness = sum(df['openness_score'])
df['A_g_max'] = df.apply(lambda row: row['openness_score'] / max_openness, axis=1)
df['A_g_sum'] = df.apply(lambda row: row['openness_score'] / sum_openness, axis=1)
df.to_csv(OUT_NAME+'.csv', float_format='%.6f')

# Create histograms
if PLOT_OPENNESS_DISTS:
    plt.figure(figsize=(8, 6))
    plt.hist(openness_hist_data, bins='auto')
    plt.xlabel('Openness score (log10)')
    plt.ylabel('# of genes')
    plt.figtext(0.99, 0.01, str(len(openness_hist_data)) + ' transcripts with peaks\n' + str(no_peaks)+' transcripts without peaks', horizontalalignment='right')
    plt.savefig(OUT_NAME+'.png')

    plt.clf()
    plt.hist(openness_gene_hist_data, bins='auto')
    plt.xlabel('Openness score (log10)')
    plt.ylabel('# of genes')
    plt.figtext(0.99, 0.01, str(len(openness_gene_hist_data))+' transcripts with peaks in gene body\n'+str(len(openness_hist_data))+' transcripts with peaks', horizontalalignment='right')
    plt.savefig(OUT_NAME+'_gene.png')

    plt.clf()
    plt.hist(openness_promoter_hist_data, bins='auto')
    plt.title('Length of promoter region = ' + str(PROMOTER_BUFFER))
    plt.xlabel('Openness score (log10)')
    plt.ylabel('# of genes')
    plt.figtext(0.99, 0.01, str(len(openness_promoter_hist_data))+' transcripts with peaks in promoter\n'+str(len(openness_hist_data))+' transcripts with peaks', horizontalalignment='right')
    plt.savefig(OUT_NAME+'_promoter.png')

if PLOT_PVAL_DIST:
    plt.clf()
    plt.hist(p_hist_data, bins='auto')
    plt.xlabel('p-value')
    plt.ylabel('# of genes')
    plt.savefig(OUT_NAME+'_p.png')

if PLOT_ACCESSIBILITY_DISTS:
    plt.clf()
    plt.hist([math.log10(p) for p in list(df['A_g_max']) if p != 0], bins='auto')
    plt.xlabel('Probability of accessibility (log10)')
    plt.ylabel('# of genes')
    plt.savefig(OUT_NAME+'_aMax.png')

    plt.clf()
    plt.hist([math.log10(p) for p in list(df['A_g_sum']) if p != 0], bins='auto')
    plt.xlabel('Probability of accessibility (log10)')
    plt.ylabel('# of genes')
    plt.savefig(OUT_NAME+'_aSum.png')

# Volcano plot bewteen 
if VOLCANO_PLOT:
    visuz.gene_exp.volcano(df=volcano_df, lfc='log2fc', pv='p-value', figname=OUT_NAME+'_volcano')