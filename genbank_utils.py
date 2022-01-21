#genbank_utils.py

from Bio import SeqIO

from override import GENE_NAME_OVERRIDE, GENE_PRODUCT_OVERRIDE


# feature tuple indices
LEFT_IDX = 0
RIGHT_IDX = 1
STRAND_IDX = 2
LOCUS_IDX = 3
GENE_IDX = 4
TYPE_IDX = 5


def get_genome_from_genbank(gb_file):
    '''
    Given a genbank file, return its underlying DNA sequence
    '''
    seq_record = SeqIO.parse(gb_file, "genbank").__next__()
    return seq_record.seq

def get_genome_fwd_rev_len(gb_file):
    '''
    Return the genome strings for the fwd and reverse strands
    and the length
    '''
    genome = get_genome_from_genbank(gb_file)
    genome_fwd = str(genome)
    genome_rev = str(genome.reverse_complement())
    genome_len = len(genome_fwd)

    return genome_fwd, genome_rev, genome_len


def get_feature_tuples_from_genbank(gb_file):
    '''
    Given a genbank file, parse out all of it's features into a 5-tuple 
    of (start_coord, end_coord,locus_tag,gene_symbol,type).
    '''
    # Use BioPython genbank parser
    seq_record = SeqIO.parse(gb_file, "genbank").__next__()
    
    feat_list = []
    # Loop over the genome file, get the CDS features on each of the strands
    for feature in seq_record.features:
        if 'locus_tag' in feature.qualifiers:
            # get  locus tag 
            lt = feature.qualifiers['locus_tag'][0]
            # get the gene symbol if available, otherwise leave blank
            g = "" if 'gene' not in feature.qualifiers else feature.qualifiers['gene'][0]
            
            feat_list.append((feature.location.start.position,
                             feature.location.end.position,
                             feature.strand,
                             lt,
                             g,
                             feature.type))
            
    return feat_list

def get_gbfeats2exclude():
    '''
    Open the config file specifying which kinds of genbank features
    to exclude when extracting pos and neg feats (these are typically
    features that a user may find are just redundant entries. For example,
    in the M. buryatense genbank file, there are duplicate entries labeled
    as "CDS" and "gene" that have the same coordinate/locus tag. I decided
    to go with "CDS" and exclude the "gene" version of each feature. 
    
    A user can amend this configuration in config/gbfeats2exclude.txt"
    '''
    with open('config/gbfeats2exclude.txt','r') as f:
        lines = [x.strip() for x in f.readlines()]
    return lines


def get_pos_neg_features(gb_file):
    '''
    Given a genbank file, load its features then sort them into
    lists of positive strand feats and negative strand feats
    '''
    feats = get_feature_tuples_from_genbank(gb_file)
    gbfeats2exclude = get_gbfeats2exclude()

    feats_filt = [x for x in feats if (x[TYPE_IDX] not in gbfeats2exclude)]
        
    # separate pos and neg lists and sort by gene start coordinate
    pos_feats = [x for x in feats_filt if x[STRAND_IDX]== 1]
    pos_feats.sort(key=lambda x: x[LEFT_IDX])
    
    neg_feats = [x for x in feats_filt if x[STRAND_IDX]==-1]
    neg_feats.sort(key=lambda x: x[RIGHT_IDX])
    
    # make sure we still keep all features
    assert(len(feats_filt) == len(pos_feats)+len(neg_feats))

    return pos_feats, neg_feats


def get_feats_and_genome(gb_file):
    '''
    Return feature tuples and underlying genome from genbank file
    '''
    feats = get_feature_tuples_from_genbank(gb_file)
    # drop duplicated 'gene' wrapper
    feats_filt = [x for x in feats if x[TYPE_IDX] != 'gene']

    genome = get_genome_from_genbank(gb_file)
    return feats_filt, genome


def get_override_gene(locus_tag,cur_gene):
    '''
    Given a locus tag, return an overridden gene
    '''
    return GENE_NAME_OVERRIDE[locus_tag] if locus_tag in GENE_NAME_OVERRIDE else cur_gene

def get_override_product(locus_tag,cur_prod):
    '''
    Given a locus tag, return an overridden gene
    '''
    return GENE_PRODUCT_OVERRIDE[locus_tag] if locus_tag in GENE_PRODUCT_OVERRIDE else cur_prod


def get_feat2meta_dict(genbank_path):
    '''
    Given a genbank file, parse it and return a dictionary of locus and the 
    gene, product and type fields
    '''
    seq_record = SeqIO.parse(genbank_path, "genbank").__next__()
    feat_list = []
    # Loop over the genome file, get the features on each of the strands
    for feature in seq_record.features:
        if feature.type != 'gene': # exclude 'gene' wrapper type
            if 'locus_tag' in feature.qualifiers: # exclude features without a locus tag
                # get  locus tag, feature name and product
                lt = feature.qualifiers['locus_tag'][0]
                g = "" if 'gene' not in feature.qualifiers else feature.qualifiers['gene'][0]
                prod = "" if 'product' not in feature.qualifiers else feature.qualifiers['product'][0]
                t = feature.type
                strand = feature.strand

                # overrides
                g = get_override_gene(lt,g)
                prod = get_override_product(lt,prod)

                metadata = {
                    'gene_symbol':g,
                    'product':prod,
                    'type':t,
                    'strand':strand
                }

                feat_list.append((lt,metadata))

    return dict(feat_list)


def get_pos_neg_relative_feature_coords(gb_file, genome_len):
    '''
    Given a genbank file, load its features and return 
    a simple tuple of it's (start, end, locus_tag). For
    feats on the negative strand, give it's coords relative
    to the genome reverse complement.
    '''

    pos_feats, neg_feats = get_pos_neg_features(gb_file)

    pos_feat_coords = [(x[LEFT_IDX], 
                        x[RIGHT_IDX],
                        x[LOCUS_IDX]) for x in pos_feats]
    
    # subtract from genome len to get position on rev compelement
    # also reverse the list
    neg_feat_coords = [(genome_len - x[RIGHT_IDX], 
                        genome_len - x[LEFT_IDX],
                        x[LOCUS_IDX]) for x in neg_feats][::-1]

    return pos_feat_coords, neg_feat_coords








