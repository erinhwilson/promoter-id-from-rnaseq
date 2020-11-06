from collections import Counter
import logomaker
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import random
random.seed(100)
import seaborn as sns
import swifter


import Bio
from Bio import motifs
from Bio.Seq import Seq

# +-----------------------------------+
# | Functions for Viewing motif logos |
# +-----------------------------------+

def load_promoter_seqs(filename):
    '''
    Load fasta file of promoters into ID, desc, and seq. It expects
    each fasta header to be divided by "|" with in the format:
    LOCUS_TAG|GENE_SYMBOL|PRODUCT
    '''
    proms = []
    with open(filename,'r') as f:
        for line in f:
            if line.startswith(">"):
                full_header = line.strip()[1:].strip()
                locus_tag = full_header.split('|')[0]
            else:
                seq = line.strip().upper()
                proms.append((locus_tag,full_header,seq))
                
    return proms


def view_motif(m1,m2):
    '''
    Given two Motif objects from the BioPython motifs module,
    create a sequence logo from their PWMs
    '''
    df1 = pd.DataFrame(m1.pwm, columns=['A','C','G','T'])
    df2 = pd.DataFrame(m2.pwm, columns=['A','C','G','T'])

    # initialize 2 panel figure
    fig,(ax1,ax2) = plt.subplots(1,2,sharey=True,figsize=[10,2])

    # block 1
    logo1 = logomaker.transform_matrix(df1,from_type='probability',to_type='information')
    crp_logo1 = logomaker.Logo(logo1,ax=ax1)

    # block 2
    logo2 = logomaker.transform_matrix(df2,from_type='probability',to_type='information')
    crp_logo2 = logomaker.Logo(logo2,ax=ax2,)

    # labels
    df1_title = f"Motif Block 1"
    df2_title = f"Motif Block 2"
    ax1.set_title(df1_title)
    ax2.set_title(df2_title)
    ax1.set_xticks([])
    ax2.set_xticks([])

    plt.ylim(0,2)

    plt.show()


def build_2Bmotif_from_selection_file(filename, randomize=False,verbose=False):
    '''
    Given a SELECTION.fa, build a 2 block consensus motif from the 
    first 6 and last 6 bases of each input sequence
    '''
    proms = load_promoter_seqs(filename)

    # collect first 6 and last 6 bases from each predicted promoter
    block1_instances = [Seq(x[:6]) for (_,_,x) in proms]
    block2_instances = [Seq(x[-6:]) for (_,_,x) in proms]
    # if we are randomizing the sequences to make a random motif
    if randomize:
        # establish a new random ordering of indicies (some may be repeated or left out)
        new_order = [random.randint(0,5) for x in range(6)]
        # remap letters to different bases. This tries to maintain a "strong signal" 
        # exists (just shuffling all the letters to different positions independently
        # would make the information content just go flat). Instead, this just changes 
        # which letters have stronger or weaker signal, but the signal strength of the 
        # original consensus remains.... I think.
        remap = {
            'A':'T',
            'C':'A',
            'G':'C',
            'T':'G',
        }
        # collect randomized versions of the hexamers
        rand_block1_instances = []
        rand_block2_instances = []
        for i in range(len(block1_instances)):
            seq1 = str(block1_instances[i])
            seq1 = ''.join([remap[x] for x in seq1])
            seq1 = ''.join([seq1[x] for x in new_order])
            rand_block1_instances.append(Seq(seq1))

            seq2 = str(block2_instances[i])
            seq2 = ''.join([remap[x] for x in seq2])
            seq2 = ''.join([seq2[x] for x in new_order])
            rand_block2_instances.append(Seq(seq2))

        block1_instances = rand_block1_instances
        block2_instances = rand_block2_instances


    # create BioPython motif objects
    m1 = motifs.create(block1_instances)
    m2 = motifs.create(block2_instances)

    # add pseudocount for proper pssm
    m1.pseudocounts = 0.5
    m2.pseudocounts = 0.5
    
    if verbose:
        # Display Motif matrix info and consensus
        print("PWM")
        print(m1.pwm)
        print(m2.pwm)
        print("\nPSSM")
        print(m1.pssm)
        print(m2.pssm)
        print("\nConsensus")
        print(m1.consensus)
        print(m2.consensus)

    # view the consensus from this file
    view_motif(m1,m2)
    
    return proms, m1, m2


# +-----------------------------------------------------+
# | Functions for scanning the genome for motif matches |
# +-----------------------------------------------------+

def build_feature_distance_index(feat_coords, genome_len):
    '''
    Given a list of feature coords on a particular strand, go through the genome and for
    each position, record the distance to the next feature start position. If the position
    is inside a feature, mark as -1.

    feat_coords should be alist of tuples: (feat_start_coord, feat_end_coord, locus_id)
    '''
    start_idx = 0
    end_idx = 1
    loc_idx = 2

    feat_idx = 0 # current feature index
    cur_feat = feat_coords[feat_idx]
    
    # keep track of distance to next feat as well as 
    # which feat we're that distance from
    dist_array = np.zeros(genome_len,dtype='U15')
    nearest_feat_array = np.zeros(genome_len,dtype='U15')
    
    # for each position in the genome
    for genome_idx in range(genome_len):
        
        # if we're before the next feature
        if genome_idx < cur_feat[start_idx]:
            # record distance from current loc to feature start
            dist_array[genome_idx] = cur_feat[start_idx] - genome_idx
        
        # if we're in the feature, mark as -1
        elif genome_idx >= cur_feat[start_idx] and genome_idx <= cur_feat[end_idx]:
            dist_array[genome_idx] = -1
        
        # if we've passed through the current feature and are on to the next, 
        # update the current feature
        elif genome_idx > cur_feat[end_idx]:
            feat_idx += 1
            # if we've run out of features at the end of the genome
            if feat_idx >= len(feat_coords):
                # build the final feat around the circle of the genome
                cur_feat = (genome_len+feat_coords[0][start_idx], genome_len+feat_coords[0][end_idx], feat_coords[0][loc_idx])
            else:
                cur_feat = feat_coords[feat_idx]

            dist_array[genome_idx] = cur_feat[start_idx] - genome_idx
        
        # mark nearest gene
        nearest_feat_array[genome_idx] = cur_feat[loc_idx] # id of feature
            
    return dist_array, nearest_feat_array


def build_genome_position_category_df(pos_dist_arr, neg_dist_arr):
    '''
    Given an array of genome positions and their distance to the 
    next nearest feature, build a dictionary counting the number 
    of positions in each type of category:
        * Inside a gene
        * within 100bp of a gene start
        * between 300bp and 100bp of a gene start
        * intergenic, beyond 300 bp from gene start

    This df is later used to normalize pssm match counts (by the total
    number of positions in the genome that are in each category)
    '''
    def get_category(num):
        '''
        Given a distance, return it's category  
        '''
        if num == -1:
            return "in gene"
        elif num <= 100:
            return "<100 to ATG"
        elif num <=300:
            return "100:300 to ATG"
        else:
            return "intergenic"

    # count number of positions in each genome category
    pos_cat_dict = Counter([get_category(int(x)) for x in pos_dist_arr])
    neg_cat_dict = Counter([get_category(int(x)) for x in neg_dist_arr])
    
    # convert to df
    pos_cat_df = pd.DataFrame.from_dict(pos_cat_dict,orient='index').rename(columns={0:'pos_count'})
    neg_cat_df = pd.DataFrame.from_dict(neg_cat_dict,orient='index').rename(columns={0:'neg_count'})

    # combine pos and negative dfs, calculate total
    cat_df = pd.concat((pos_cat_df,neg_cat_df),axis=1).reset_index().rename(columns={'index':'cat'})
    cat_df['total'] = cat_df.apply(lambda row: row['pos_count'] + row['neg_count'],axis=1)
    
    # add row combining results from both <100 and <300 of ATG
    all_300_pos = pos_cat_dict['100:300 to ATG'] + pos_cat_dict['<100 to ATG']
    all_300_neg = neg_cat_dict['100:300 to ATG'] + neg_cat_dict['<100 to ATG']
    all_300_total = all_300_pos + all_300_neg
    row = pd.DataFrame([['<300 to ATG',all_300_pos, all_300_neg, all_300_total]],
                       columns=['cat','pos_count','neg_count','total'])
    cat_df = pd.concat((cat_df,row),ignore_index=True)
    
    return cat_df

# go through all the pssm matches found and determine they're intergenic category
def get_intergenic_category(row,dist_array):
    '''
    Given a pssm match row and a dist array, return the type of 
    intergenic category it is based on distance to the feature
    '''
    # get the ending position of the motif (2 hexamers + spacer -1 for indexing)
    pssm_start_pos = row['pos']
    pssm_end_pos = pssm_start_pos + 12 + row['spacer'] - 1
    
    pos_dist = int(dist_array[pssm_end_pos])
    
    if pos_dist == -1:
        cat = "in gene"
    elif pos_dist <= 100:
        cat = "<100 to ATG"
    elif pos_dist <=300:
        cat = "100:300 to ATG"
    else:
        cat = "intergenic"
        
    return cat
    
    
def add_intergenic_category_column(row, 
                                   pos_dist_array, 
                                   neg_dist_array):
    '''
    Given a pssm match row, get it's intergenic category. Take into account fwd vs rev
    '''
    genome_version = row['seq_id']
    
    if genome_version == 'genome_fwd':
        return get_intergenic_category(row,pos_dist_array)
    elif genome_version == 'genome_rev':
        return get_intergenic_category(row,neg_dist_array)
    else: 
        raise ValueError(f"Unknown genome direction {genome_version}. This function is for whole genome motif searches (genome_fwd or genome_rev)")

def add_nearest_feat_column(row,pos_nearest_feat_array,neg_nearest_feat_array):
    '''
    Given a pssm match row, get it's nearest feature. Take into account fwd vs rev
    '''
    genome_version = row['seq_id']
    
    if genome_version == 'genome_fwd':
        return pos_nearest_feat_array[row['pos']]
    elif row['seq_id'] == 'genome_rev':
        return neg_nearest_feat_array[row['pos']]
    else: 
        raise ValueError(f"Unknown genome direction {genome_version}. This function is for whole genome motif searches (genome_fwd or genome_rev)")

        
genome_categories = ['in gene','intergenic','100:300 to ATG','<100 to ATG']
genome_cat_order = ['in gene','intergenic','<300 to ATG','<100 to ATG']
def genome_category_violin(df,threshold=0):
    fig = plt.figure(figsize=(10,10))
    sns.violinplot(data=df[df['score']>=threshold],
                   x='motif_loc',
                   y='score',
                   order = genome_cat_order
                  )
    plt.xlabel("Genome category")
    plt.ylabel("PSSM match score")
    plt.title(f"Distribution of PSSM match scores (>={threshold})")
    
    plt.show()
    
def genome_category_swarm(df,threshold=0):
    fig = plt.figure(figsize=(20,10))
    sns.swarmplot(data=df[df['score']>threshold], 
                  x='motif_loc',
                  y='score',
                  order = genome_cat_order,
                  hue='spacer')
    plt.xlabel("Genome category")
    plt.ylabel("PSSM match score")
    plt.title(f"Distribution of PSSM match scores (>={threshold})")
    plt.show()
    
def genome_category_normed_bar_v(df,threshold=0,sci=False):
    plt.figure(figsize=(4,7))
    sns.barplot(data=df, x='cat',
                y='match_perc',
                order = genome_cat_order
               )
    plt.yticks(fontsize=14)
    plt.xticks(rotation=90,fontsize=14)
    if sci:
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    plt.xlabel("Genome category",fontsize=14)
    plt.ylabel("Fraction of genome positions with PSSM match",fontsize=14)
    plt.title(f"Enrichment of PSSM matches \nin each genome category \n(PSSM log odds >{threshold:0.2f})",fontsize=20)
    plt.show()

def genome_category_normed_bar_h(df,threshold=0,sci=False):
    plt.figure(figsize=(8,4))
    sns.barplot(data=df, y='cat',
                x='match_perc',
                order = genome_cat_order
               )
    plt.yticks(fontsize=14)
    plt.xticks(fontsize=14)
    if sci:
        plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
    plt.ylabel("Genome category",fontsize=14)
    plt.xlabel("Fraction of genome positions with PSSM match",fontsize=14)
    plt.title(f"Enrichment of PSSM matches \nin each genome category \n(PSSM log odds >{threshold:0.2f})",fontsize=20)
    plt.show()
    
def add_genome_category_to_pssm_matches(df,
                                        pos_dist_array,
                                        pos_nearest_feat_array,
                                        neg_dist_array,
                                        neg_nearest_feat_array,
                                        make_swarm_violin=False):

    '''
    Given a df of pssm motif match scores across the entire genome, 
    Add the genome category and teh nearest feature to the dataframe 
    of pssm matches 
    '''
    # add genome pos category to each match in the df
    print("Adding categories to pssm matches...")
    df['motif_loc'] = df.swifter.apply(lambda row: add_intergenic_category_column(row, pos_dist_array, neg_dist_array),axis=1)
    df['nearest_feat'] = df.swifter.apply(lambda row: add_nearest_feat_column(row,pos_nearest_feat_array,neg_nearest_feat_array),axis=1)

    if make_swarm_violin:
        # make a violinplot
        print("Plotting violin and swarm...")
        genome_category_violin(df)
        genome_category_violin(df,threshold=thresh1)

        # make a swarm plot above threshold
        genome_category_swarm(df,threshold=thresh2)

    return df 

def analyze_motif_matches_across_genome(df,
                                        baseline_cat_df):
    '''
    Given a df of pssm motif match scores across the entire genome, 
    analyze the frequency of these matches in various positions relative
    to genes. Are the matches mostly in genes? in promoter regions? 
    Intergenic but far away from promoter regions? 
    '''
    

    # now normalize the match counts by the number of positions in the genome that
    # fall into each category
    print("Normalizing pssm match counts")

    # chekc if df is empty
    # if df.empty:
    #     pssm_cat_match_counts = dict([(x,0) for x in genome_categories])
    # else:
        
    # make a dictionary of each genome category and the number of matches
    pssm_cat_match_counts = dict(df['motif_loc'].value_counts())
    # catch any empty categories
    for cat in genome_categories:
        if cat not in pssm_cat_match_counts:
            pssm_cat_match_counts[cat] = 0
    
    # combine <100 and 100:300 to be all <300
    # c100 = 0 if '<100 to ATG' not in pssm_cat_match_counts else pssm_cat_match_counts['<100 to ATG']
    # c100_300 = 0 if '100:300 to ATG' not in pssm_cat_match_counts else pssm_cat_match_counts['100:300 to ATG']
    
    pssm_cat_match_counts['<300 to ATG'] = pssm_cat_match_counts['<100 to ATG'] + \
                                           pssm_cat_match_counts['100:300 to ATG']

    
    cat_df = baseline_cat_df.copy(deep=True)
    
    # add pssm match count to this df
    cat_df['pssm_match_count'] = cat_df['cat'].apply(lambda x: pssm_cat_match_counts[x])
    # divide match count by total counts of each category
    cat_df['match_perc'] = cat_df.apply(lambda row: row['pssm_match_count']/row['total'], axis=1)

    return cat_df


def make_spaced_motif(m1,m2,spacer):
    '''
    Given 2 pssm blocks, build a full pssm connecting them
    by the specified spacer distance
    '''
    full_motif = {}
    
    for base in ['A','C','G','T']:
        full_motif[base] = m1.pssm[base] + \
                           [0.0 for x in range(spacer)] + \
                           m2.pssm[base]
    full_motif = Bio.motifs.matrix.PositionSpecificScoringMatrix(m1.alphabet,full_motif)
    return full_motif


def build_dict_of_motifs_to_try(m1,
                                m2,
                                spacers=[15,16,17,18]
                               ):
    '''
    Given a list of spacers to try, build a dictionary of motif PSSM objects
    '''
    motifs_to_try = {}

    for sp in spacers:
        motifs_to_try[sp] = make_spaced_motif(m1,m2,sp)

    return motifs_to_try


def find_and_score_motifs_in_seqs(motif_dict,seqs,truth_dict):
    '''
    Given a dictionary of motifs (key:spacer, value:biopython motif pssm),
     a list of promoter sequences in which to search for the motifs, and a 
     dictionary of our best ground truth guess for the the best actual motif,
     use BioPython motif search tools to find and score motifs in the sequences.
     Make a note if the found sequence exactly matches the ground truth
    '''
    
    # for each promoter
    rows = []
    for seq_id,seq_name,seq in seqs:
        # for each spacer version of the motif
        for spacer in motif_dict:
            # get the pssm object from the dict
            pssm = motif_dict[spacer]
            # for each result from searching the seq for the pssm
            try:
                search_res = pssm.search(seq,both=False)
                
                for pos,score in search_res:
                    extracted_seq = seq[pos:pos+len(pssm['A'])]
                    # does this match the ground truth from BioP?
                    if seq_id in truth_dict:
                        match_best = truth_dict[seq_id] == extracted_seq
                    # if the seq is not in the truth dict, then we don't have a guess for it. default false
                    else:
                        match_best = False

                    rows.append([seq_id,seq_name,score,pos,len(seq),spacer,extracted_seq,match_best])
            
            except ValueError as e:
                print(f"{seq_id} and {spacer}-spacer PSSM search error: {e}")
                print("skipping...")            
            
        
    df = pd.DataFrame(rows,columns=['seq_id','seq_name','score','pos','seq_len','spacer','full_seq','match_best?'])
    
    return df


def score_predictions_to_motif(motif_blocks, m1, m2):
    '''
    Given a 2 block consensus motif (m1 and m2), as the predicted sequences 
    that comprise the consensus motifs, score each sequence against the 
    consensus pssm.
    '''

    # for each promoter prediction, score it's hexamers relative to the consensus
    hex_score_data = []
    for loc,desc,seq in motif_blocks:
        hex1 = seq[:6]
        hex2 = seq[-6:]
        hex1_score = m1.pssm.calculate(hex1)
        hex2_score = m2.pssm.calculate(hex2)
        total_score = hex1_score + hex2_score

        row = [loc,desc,seq,hex1,hex1_score,hex2,hex2_score,total_score]
        hex_score_data.append(row)
        
    hex_score_df = pd.DataFrame(hex_score_data,columns=['locus_tag','desc','motif_block','hex1','hex1_score','hex2','hex2_score','total_score'])

    return hex_score_df







