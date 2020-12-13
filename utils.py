import pandas as pd
import copy

"""
Helper functions and utilities that may be useful for future analyses
- FBgn -> gene symbol conversion
- Operations for vector representation of base pair sequences
- Ratcliff-Obershelp sequence similarity
"""

# Convert an FBgn to a gene symbol
conversions = pd.read_csv('fbgn_annotation_ID_fb_2020_04.tsv', sep='\t', header=0)
def FBgn_to_genesymbol(FBgn):
    try:
        genesymbol = conversions[conversions['primary_FBgn#'] == FBgn]['gene_symbol'].iloc[0]
    except:
        try:
            genesymbol = conversions[conversions['secondary_FBgn#(s)'] == FBgn]['gene_symbol'].iloc[0]
        except:
            genesymbol = FBgn
    return genesymbol

"""
BP Vector Operations
 Vector representation: a one-hot vector, where each group of 4 elements represents one
 base pair position in a human-readable sequence
 e.g. 'ATTG' -> 1000 0100 0100 0001
"""
A = [1, 0, 0, 0]
T = [0, 1, 0, 0]
C = [0, 0, 1, 0]
G = [0, 0, 0, 1]

def seq_to_vec(seq):
    vec = []
    for base in seq:
        if base == 'A':
            vec += A
        if base == 'T':
            vec += T
        if base == 'C':
            vec += C
        if base == 'G':
            vec += G
    return vec

def vec_to_seq(vec):
    seq = ''
    for i, val in enumerate(vec):
        if val == 2:
            seq += '-'
        if val != 1:
            continue
        if i%4 == 0:
            seq += 'A'
        elif i%4 == 1:
            seq += 'T'
        elif i%4 == 2:
            seq += 'C'
        elif i%4 == 3:
            seq += 'G'
    return seq

# Converts a vector of floats to a base pair sequence representation by choosing
# the maximum value out of each 4-element group.
# If the maximum does not meet the threshold, that position will be represented by a -
def floatvec_to_seq(floatvec, threshold=0):
    quads = [floatvec[i:i + 4] for i in range(0, len(floatvec), 4)]
    vec = []
    for quad in quads:
        quad = [abs(val) for val in quad]
        one = list(quad).index(max(quad))
        one_hot = [0, 0, 0, 0]
        one_hot[one] = 1 if max(quad) >= threshold else 2
        vec += one_hot
    return vec_to_seq(vec)

# Adds fuzziness to a vector representation by spreading each occurance of a 1
# e.g. 1000 0100 0100 0001 -> 1,.25,0,0  .25,1,.25,0 .25,1,.25,0  0,0,.25,1
#                       OR -> 1100 1110 1110 0011 (if course=True) 
def spread_vec(vec, course=False):
    spread_factor = 1 if course else 0.25
    out = copy.deepcopy(vec)
    for i, bp in enumerate(vec):
        if bp == 0:
            continue
        if i-1 >= 0 and out[i-1] < 1:
            out[i-1] += bp*spread_factor
        if i+1 < len(vec)-1 and out[i+1] < 1:
            out[i+1] += bp*spread_factor
    return out

"""
Ratcliff-Obershalp Similarity / Gestalt Pattern Matching
Measures sequence similarity by comparing the length of the longest common substrings
to the length of the two strings
"""
# https://www.geeksforgeeks.org/python-program-for-longest-common-subsequence/#:~:text=Python%20Program%20for%20Longest%20Common%20Subsequence.%201%20filter_none.,1%29%3A%20if%20i%20%3D%3D%200%20...%20More%20items
def longest_common_seq(X, Y):
    # find the length of the strings 
    m = len(X) 
    n = len(Y) 
  
    # declaring the array for storing the dp values 
    L = [[None]*(n + 1) for i in range(m + 1)] 
  
    """Following steps build L[m + 1][n + 1] in bottom up fashion 
    Note: L[i][j] contains length of LCS of X[0..i-1] 
    and Y[0..j-1]"""
    for i in range(m + 1): 
        for j in range(n + 1): 
            if i == 0 or j == 0 : 
                L[i][j] = (0, '')
                #print(L[i][j], 1)
            elif (X[i-1] == Y[j-1]) and (L[i-1][j-1][1]+X[i-1] in X and L[i-1][j-1][1]+X[i-1] in Y): 
                L[i][j] = (L[i-1][j-1][0]+1, L[i-1][j-1][1]+X[i-1])
                #print(L[i][j], 2, X[i-1], Y[j-1])
            else: 
                if L[i-1][j][0] > L[i][j-1][0]:
                    L[i][j] = L[i-1][j]
                else:
                    L[i][j] = L[i][j-1]
                #print(L[i][j], 3, L[i][j][1] in X and L[i][j][1] in Y)
  
    # L[m][n] contains the length of LCS of X[0..n-1] & Y[0..m-1] 
    return L[m][n] 

# Recursion step
def common_substr_length(p1, p2):
    if len(p1) == 0 or len(p2) == 0:
        return 0
    if len(p1) == 1 and len(p2) == 1:
        return 1 if p1 == p2 else 0
    
    lcs_len, lcs = longest_common_seq(p1, p2)
    
    p1_start = p1.index(lcs)
    p2_start = p2.index(lcs)
    
    p1_left = p1[:p1_start]
    p2_left = p2[:p2_start]
    p1_right = p1[p1_start+lcs_len:]
    p2_right = p2[p2_start+lcs_len:]
    
    ll_similarity = common_substr_length(p1_left, p2_left)
    rr_similarity = common_substr_length(p1_right, p2_left)
    lr_similarity = common_substr_length(p1_left, p2_right)
    rl_similarity = common_substr_length(p1_right, p2_left)
    
    return lcs_len + (ll_similarity+lr_similarity)/2 + (rr_similarity + rl_similarity)/2 
    
# Ratcliff-Obershelp similarity
def sequence_similarity(p1, p2):
    similar_chars = common_substr_length(p1, p2)
    return (2*similar_chars)/(len(p1)+len(p2))