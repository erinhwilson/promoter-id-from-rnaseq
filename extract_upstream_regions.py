# extract_upstream_regions.py
# script that accepts a 2 column file of locus
# ids and a boolean flag of whether than locus
# might be in an operon. It returns a fasta file
# from 1 base before the ATG through -300bp
# upstream for all genes NOT possibly in an operon.


import argparse
from Bio import SeqIO
import os
import pandas as pd

import genbank_utils as gu
import get_top_gene_set as gtgs 

# genbank feature tuple indices
LEFT_IDX = 0
RIGHT_IDX = 1
STRAND_IDX = 2
LOCUS_IDX = 3
GENE_IDX = 4
TYPE_IDX = 5

def load_loci(filename):
    df = pd.read_csv(filename,sep='\t')

    # return list of loci not flagged as potentially in an operon
    loci = list(df[df['op?']==False]['locus_tag'].values)
    return loci 


def get_all_upstream_regions(gb_file, 
                             window_size, 
                             min_dist,# default=20 in argparse
                             no_truncate,# default=False in argparse
                             avoid_rbs,# default=None in argpase
                             pre_window_extension,#default=0 in argparse
                             verbose=False,
                             ):
    '''
    Given a genbank file, parse out all its features into 6-tuples
    of (LEFT, RIGHT, STRAND, LOCUS, GENE, TYPE). Then use the sequence
    in the genbank file to extract a window_size number of base pairs
    upstream of each features start coordinate.

    Truncate is on by default and will stop slicing if we run into 
    the coding sequence of another feature (so only extract the 
    intergenic sequence). 

    However some features are too close together (and possibly even 
    overlapping). While most of these should have been ignored in 
    the operon estimation step of get_top_gene_sets.py, some genes
    are oriented divergently and thus still not in an operon. If
    the distance between two loci is < min_dist, we extract min dist
    anyways, even if it runs into another locus annotation.
    '''
    # load genbank features and genome
    feats,genome = gu.get_feats_and_genome(gb_file)

    # dictionary to collect feature upstream seqs
    upstream_regions = {}

    # loop through features and get upstream for all
    for i,cur_feat in enumerate(feats):
        # +-----------------+
        # | NEGATIVE STRAND |
        # +-----------------+
        # if we're on the negative strand, go 300bp to the right, reverse compelement
        if cur_feat[STRAND_IDX] == -1:
            # get the range of the promoter region
            p_left = cur_feat[RIGHT_IDX] + 1
            p_right = p_left + window_size

            # if there's a request to extend into the start of the gene, add it here
            if pre_window_extension:
                p_left -= pre_window_extension

            # extract a slice from the genome at the proper coordinates
            # revcomp because on reverse
            seq = genome[p_left:p_right].reverse_complement()

            # if in truncate mode, check if upstream feat is too close
            if not no_truncate:
                # make sure we're not at the last feat in the list (no
                # rightward gene upstream)
                if i < len(feats) - 1:
                    # get the FOLLOING feature (because on -1 strand)
                    upstream_feat = feats[i+1]
                    # how far is the upstream feat from the current?
                    upstream_dist = upstream_feat[LEFT_IDX] - cur_feat[RIGHT_IDX]
                    # if it's closer than the window we extracted, we need
                    # to truncate    
                    if upstream_dist < window_size:
                        if verbose:
                            if upstream_dist < min_dist:
                                print("SHORT DIST!")
                                print("Cur feat:")
                                print(cur_feat)
                                print("Up feat:")
                                print(upstream_feat)
                                print(f"distance:{upstream_dist}\n")
                        # if upstream distance is too small (features are closer than
                        # min_dist), then set the upstream dist to at least min_dist
                        upstream_dist = max(upstream_dist, min_dist)

                        # determine how much of the window to truncate to avoid
                        # running into the next feature
                        trunc_dist = window_size - upstream_dist

                        # truncate from the beginning of the seq (the upstreamer part)
                        # (we've already reverse complemented so we can still take the 
                        # upstream part of the promoter seq even though this is the reverse 
                        # section)
                        seq = seq[trunc_dist:]

        # +-----------------+
        # | POSITIVE STRAND |
        # +-----------------+
        # if we're on the positive strand, go 300bp to the left
        elif cur_feat[STRAND_IDX] == 1:
            p_right = cur_feat[LEFT_IDX] - 1 
            p_left = p_right - window_size

            # if there's a request to extend into the start of the gene, add it here
            if pre_window_extension:
                p_right += pre_window_extension
            
            # extract a slice from the genome at the proper coordinates
            seq = genome[p_left:p_right]

            # if in truncate mode
            if not no_truncate:
                # make sure this isn't the very first feat in the list
                # (no leftward upstream gene)
                if i > 0:
                    # get the PREVIOUS feature (because on +1 strand)
                    upstream_feat = feats[i-1]
                    # how far is the upstream feat from the current?
                    upstream_dist = cur_feat[LEFT_IDX] - upstream_feat[RIGHT_IDX]

                    # if it's closer than the window we extracted, we need
                    # to truncate    
                    if upstream_dist < window_size:
                        if verbose:
                            if upstream_dist < min_dist:
                                print("SHORT DIST!")
                                print("Cur feat:")
                                print(cur_feat)
                                print("Up feat:")
                                print(upstream_feat)
                                print(f"distance:{upstream_dist}\n")
                        # if upstream distance is too small (features are closer than
                        # min_dist), then set the upstream dist to at least min_dist
                        upstream_dist = max(upstream_dist, min_dist)

                        # determine how much of the window to truncate to avoid
                        # running into the next feature
                        trunc_dist = window_size - upstream_dist

                        # truncate from the beginning of the seq (the upstreamer part)
                        seq = seq[trunc_dist:]

        else:
            raise ValueError(f"Unknown strand type: {cur_feat[STRAND_IDX]}. Expected 1 or -1.")


        # add the feature and its upstream seq to the dict
        # key : locus_tag, value : upstream sequence string
        if avoid_rbs:
            # if the sequence is already short and the avoid_rbs amount 
            # would reduce the seq to shorter than min_dist, only reduce
            # to the min dist. 
            bp_to_cut = min(avoid_rbs, len(seq)-min_dist)

            if verbose:
                print(f"Truncating {bp_to_cut} bases from end of seq")
            seq = seq[:-bp_to_cut]

        upstream_regions[cur_feat[LOCUS_IDX]] = seq 

    return upstream_regions  

def write_fasta_file(args, loci, upstream_regions):
    # load feature meta data
    feat2meta = gu.get_feat2meta_dict(args.gb_file)

    # construct output file name:
    # use same base name as loci file
    base = os.path.basename(args.loci_file).split('.')[0]
    # append "_trunc" if in truncation mode
    trunc_string = "" if args.no_trunc else "_trunc"
    # append rbs string if in rbs_avoidance mode
    rbs_flag = "" if not args.avoid_rbs else f"_RBSminus{args.avoid_rbs}"
    # append pre-ext string if in an extension was added
    pre_ext_flag = "" if not args.pre_window_ext else f"_preext{args.pre_window_ext}"
    # concat some relevant args
    filename = f"{base}_upstream_regions_w{args.window_size}{pre_ext_flag}{rbs_flag}_min{args.min_dist}{trunc_string}.fa"
    # path to outdir
    out_path = os.path.join(args.outdir,filename)

    # write fasta file
    with open(out_path,'w') as f:
        for loc in loci:
            # some extra metadata for output readability
            gene_symbol = feat2meta[loc]['gene_symbol']
            product = feat2meta[loc]['product']
            header = f">{loc}|{gene_symbol}|{product}"

            f.write(f"{header}\n{upstream_regions[loc]}\n")


    return out_path



# +-------------+
# | MAIN SCRIPT |
# +-------------+

def main():
    # +------+
    # | ARGS |
    # +------+
    parser = argparse.ArgumentParser(description='Extract upstream DNA regions from set of locus ids.')
    # Required args
    parser.add_argument('loci_file', help='Two-column tab-delimited file containing list of locus IDs and boolean operon flag')
    parser.add_argument('gb_file', help='Genbank file with feature annotations')
    parser.add_argument('outdir', help='Output directory where results are written')
    # Optional args
    parser.add_argument('-w', '--window_size',default=300, type=int, help='bp length of upstream region to extract')
    parser.add_argument('-p', '--pre_window_ext',default=0, type=int, help='bp length of upstream region to extract')
    parser.add_argument('-m', '--min_dist',   default=20,type=int,help='Minimum upstream distance to extract, even if features are too close.')
    parser.add_argument('-t', '--no_trunc', action='store_true',help='Turn OFF truncation mode - so always extract window_size bp, even if it overlaps with other features')
    parser.add_argument('-r', '--avoid_rbs',nargs='?',type=int,const=15, default=None,help='Turn ON RBS avoidance to truncate the end of the extracted sequence by n bases (default n=15). It will not reduce a sequence to be shorter than min_dist')
    
    args = parser.parse_args()

    # get loci
    print("Loading loci of interest...")
    loci = load_loci(args.loci_file)
    print(loci)

    # get upstream regions
    print("Getting upstream regions from genbank")
    upstream_regions = get_all_upstream_regions(
        args.gb_file,
        args.window_size, 
        args.min_dist,
        args.no_trunc,
        args.avoid_rbs,
        args.pre_window_ext)

    
    # +------+
    # | SAVE |
    # +------+
    # save to a fasta file
    print("Saving fasta of upstream regions for loci of interest...")
    out_path = write_fasta_file(args, loci, upstream_regions)    

    print(f"Output written to {out_path}")
    print("Done!")
    

if __name__ == '__main__':
    main()