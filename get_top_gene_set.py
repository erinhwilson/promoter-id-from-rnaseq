# get_top_gene_set.py
# Script that takes in a matrix of TPM expression data, 
# groups the samples into experimental conditions and takes
# the average, then saves a list loci that are in the top
# N% across all conditions listed. 

# The user must provide the data matrix, a genbank file, a value 
# for N, and an operon min distance. Other required arguments to make
# logistics smooth are an output directory and a file mapping each 
# sample name to its assigned condition. Optionally, the user can 
# provide a list of specific samples to include and a list of 
# experimental conditions. If none are provided, the program 
# will use all available. 

import argparse
import os
import pandas as pd

import genbank_utils as gu
# genbank feature tuple indices
LEFT_IDX = 0
RIGHT_IDX = 1
STRAND_IDX = 2
LOCUS_IDX = 3
GENE_IDX = 4
TYPE_IDX = 5


# +------------------------+
# | Data Loading functions |
# +------------------------+

def load_data(
        data_mat_file, 
        sample2cond_file, 
        sample_file, 
        condition_file,
        gb_file):
    '''
    Wrapper function to load data from files into relavent objects
    '''
    # load TPM data
    df = pd.read_csv(data_mat_file,sep='\t').fillna('')

    
    # load mapping from sample to condition
    with open(sample2cond_file,'r') as f:
        sample2condition = dict(x.strip().split() for x in f.readlines())

    
    # load sample to include file
    if sample_file:
        with open(sample_file,'r') as f:
            samples = list(x.strip() for x in f.readlines())
    # if none provided, just use all the samples from the sample2condition dict
    else: 
        samples = list(sample2condition.keys())

        
    # load the conditions to include file
    if condition_file:
        with open(condition_file,'r') as f:
            conditions = list(x.strip() for x in f.readlines())
    # if none provided, just use all the conditions
    else:
        conditions = list(set([sample2condition[x] for x in sample2condition]))

    # Load features from genbank
    pos_feats, neg_feats = gu.get_pos_neg_features(gb_file)


    return df, sample2condition, samples, conditions, pos_feats, neg_feats

    # CODE REVIEW: argparse error handling for opening files


def load_data_from_args(args):
    '''
    load data from argparse object
    '''
    return load_data(args.data_mat,args.sample2cond,args.samples,args.conditions,args.gb_file)


# +---------------------------------+
# | Data Transformations & analysis |
# +---------------------------------+

def transpose_tpm_data(df, samples, conditions, sample2condition):
    '''
    Given a data matrix of rows as gene loci and columns as a mix of metadata and sample ids,
    tranpose the dataframe to by samples by genes. Use the samples and conditions provided 
    by the user to filter to only the rows needed. Also add the condition as a column.
    '''
    # transpose dataframe to make experiments the rows and genes the columns
    df_T = df.set_index('locus_tag')[samples].T.reset_index().rename(columns={'index':'sample'})
    # add experimental condition column to dataframe
    df_T['exp_condition'] = df_T['sample'].apply(lambda x: sample2condition[x])
    # filter to only conditions that are requested
    df_T = df_T[df_T['exp_condition'].isin(conditions)]

    return df_T


def get_average_tpm_by_condition(df_orig,
                                 samples,
                                 conds,
                                 sample2cond,
                                 loci,
                                 add_pseudocount=True,
                                 pseudocount = 0.01
                                ):
    '''
    Given a matrix of genes by samples, first transpose it to samples 
    by genes. Then group the samples by their experimental condition.
    Finally for each gene, calculated the average and standard deviation
    '''
    # transpose data
    df_T = transpose_tpm_data(df_orig, samples, conds, sample2cond)
    
    # start with a fresh copy of the dataframe
    df = df_T.copy()
    
    # add a pseudocount to all genes' TPMs
    if add_pseudocount:
        df[loci] = df[loci] + pseudocount 
    
    # group samples by experimental condition and take the average
    df_means = df[['exp_condition']+loci]\
                    .groupby('exp_condition')\
                    .mean()\
                    .reset_index()
    
    return df_means

def flag_potential_operon_loci(pos_feats, neg_feats, min_dist):
    '''
    Given all features on the pos and neg strands, flag which ones 
    may be in an operon at the given min_dist between loci on the same strand.
    
    Return a set of flagged loci.
    '''
    
    # collect locus tags that are potential operon candidates
    operon_candidate_loci = []
    
    # +-----------------+
    # | NEGATIVE STRAND |
    # +-----------------+
    # if we're in the negative list, check the next feature to the right
    for i,cur_feat in enumerate(neg_feats):
        # if this isn't the last gene (aka no rightward operon),
        # check the min_distance window for other annotations
        if i < len(neg_feats) -1:
            # get the FOLLOWING feature (because on -1 strand)
            upstream_feat = neg_feats[i+1]

            # if the left side of the upstream gene is within min_distance, 
            # collect it
            if upstream_feat[LEFT_IDX] < (cur_feat[RIGHT_IDX] + min_dist):
                operon_candidate_loci.append(neg_feats[i][LOCUS_IDX])
    
    # +-----------------+
    # | POSITIVE STRAND |
    # +-----------------+
    # if we're in the positive list, check the next feature to the left
    for i,cur_feat in enumerate(pos_feats):
        # if this isn't the first gene (aka no leftward operon),
        # check the min_distance window for other annotations
        if i !=0:
            # get the PREVIOUS feature (because on +1 strand)
            upstream_feat = pos_feats[i-1]

            # if the right side of the upstream gene is within min_distance, 
            # collect it
            if upstream_feat[RIGHT_IDX] > (cur_feat[LEFT_IDX] - min_dist):
                operon_candidate_loci.append(pos_feats[i][LOCUS_IDX])
    
    return set(operon_candidate_loci)


def get_top_n_perc_by_condition(df, loci, n):
    '''
    Given a dataframe of experimental conditions and the average TPM of each gene,
    loop through each condition and rank the genes in descending TPM order. Keep of list
    of genes that make the top n% in every condition
    '''
    num_loci = len(loci)
    # keep track of a list of top loci for each condition
    # key = exp_condition, value = (gene,tpm)
    top_n_loci = {}
    
    # for each experimental condition row
    for i,row in df.iterrows():
        exp_cond = row.exp_condition
        # Sort all genes in the row by their tpm counts
        tpms = sorted(
                list(zip(loci, row[loci].values)), # zip locus and TPM
                key=lambda x: x[1],                # sort by TPM
                reverse=True                       # descending TPM order
        )
    
        # keep the top n% of the list
        top_n_loci[exp_cond] = tpms[:int(0.01*n*num_loci)]
        
    # keep only the genes that were in all lists via set intersection
    # drop tpm val from tuple
    top_n_loci_only = [[y[0] for y in top_n_loci[x]] for x in top_n_loci]
    # convert lists to sets
    top_n_sets = [set(x) for x in top_n_loci_only]
    # get intersection of all sets                
    top_n_loci_all_conds = set.intersection(*top_n_sets)

    return top_n_loci_all_conds#, top_n_loci, top_n_sets

# +-------------------+
# | Saving and output |
# +-------------------+

def write_loci_to_file(loci,loci_op_filt, outfile):
    '''
    Given a list of locus tags, and a list of locus_tags possbily in an operon,
    write them to a simple tab delimited file with TRUE next to possible operons
    '''
    
    with open(outfile,'w') as f:
        # set header to be "locus_tag" \tab "op?" \newline
        f.write('locus_tag\top?\n')
        
        for loc in loci:
            op = True if loc in loci_op_filt else False
            f.write(f'{loc}\t{op}\n')

# +-------------+
# | MAIN SCRIPT |
# +-------------+

def main():
    # +------+
    # | ARGS |
    # +------+
    parser = argparse.ArgumentParser(description='Get a set of top genes across conditions.')
    # Required args
    parser.add_argument('data_mat', help='tsv of TPM expression data. rows = genes, columns = samples, gene metadata')
    parser.add_argument('locus_id_col', help='The column name to use to extract unique locus IDs')
    parser.add_argument('top_n', type=int, help='Percentage threshold - get top n\% of genes across all conditions')
    parser.add_argument('operon_min_dist', type=int, help='Min distance two genes must be apart to not be considered possibly in an operon')
    parser.add_argument('sample2cond', help='Two-column tab-delimited file mapping sample names to their experimental conditions')
    parser.add_argument('gb_file', help='Genbank file with feature annotations')
    parser.add_argument('outdir', help='Output directory where results are written')
    # Optional args
    parser.add_argument('-s', '--samples', help='Filename of samples to include')
    parser.add_argument('-c', '--conditions', help='Filename of experimental conditions to include')
    
    args = parser.parse_args()


    # +-----------+
    # | LOAD DATA |
    # +-----------+
    # load expression data from files into objects
    print("Loading data...\n")
    df,\
    sample2cond,\
    samples,\
    conds,\
    pos_feats,\
    neg_feats = load_data_from_args(args)
    n = args.top_n
    min_dist - args.operon_min_dist

    # +-------------------+
    # | OPERON ESTIMATION |
    # +-------------------+
    # flag loci potentially in operons (set of locus ids)
    print("Estimating operons... \n")
    maybe_operon_loci = flag_potential_operon_loci(pos_feats, 
                                                   neg_feats, 
                                                   min_dist)


    # +-----------+
    # | TPM MEANS |
    # +-----------+
    print("Calculating mean TPMs for each condition...\n")
    # use this column to get a full list of all genes for which expression was measured
    locus_ids = list(df[args.locus_id_col].values)

    # group by condition and calculate TPM means
    df_means = get_average_tpm_by_condition(df,
                                            samples,
                                            conds,
                                            sample2cond,
                                            locus_ids)

    
    # +--------------------+
    # | IDENTIFY TOP GENES |
    # +--------------------+
    # get top locus ids in all conditions
    print("Getting top genes across conditions...\n")
    top_locs = get_top_n_perc_by_condition(df_means, locus_ids, n)
    
    # filter out loci maybe in operons
    top_locs_op_filter_out = [x for x in top_locs if x in maybe_operon_loci]
    print(f'Top {n}%: {len(top_locs)} ({len(top_locs) - len(top_locs_op_filter_out)} filtered)\n')


    # +------+
    # | SAVE |
    # +------+
    # save to a two-column txt file
    outf = os.path.join(args.outdir,f'loci_in_top_{n}perc.txt')
    write_loci_to_file(top_locs,top_locs_op_filter_out, outf)   

    print(f"Output written to {outf}")
    print("Done!")
    

if __name__ == '__main__':
    main()



