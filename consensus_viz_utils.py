import logomaker
import matplotlib.pyplot as plt
import pandas as pd

import Bio
from Bio import motifs
from Bio.Seq import Seq

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
                seq = line.strip()
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


def build_2Bmotif_from_selection_file(filename, verbose=False):
    '''
    Given a SELECTION.fa, build a 2 block consensus motif from the 
    first 6 and last 6 bases of each input sequence
    '''
    proms = load_promoter_seqs(filename)

    # collect first 6 and last 6 bases from each predicted promoter
    block1_instances = [Seq(x[:6]) for (_,_,x) in proms]
    block2_instances = [Seq(x[-6:]) for (_,_,x) in proms]

    # creat BioPython motif objects
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