# predict_promoter_signal.py

# given a fasta file of upstream regions for a set of 
# genes, run these seqs through bioprospector to search
# for promoter motifs. By default, we search for motifs
# of the structure: hexamer , 15-18bp spacer , hexamer
# and bioprospector is run 20 times then summarized to 
# pick the "most popular" signal location for each 
# upstream sequence in the fasta.

import argparse
import os
import pandas as pd
from subprocess import Popen, PIPE
import time

import bioprospector_utils as bu 


def parse_bioprospector_config(filename):
    '''
    Parse the arguments out of the bioprospector config file
    '''
    with open(filename,'r') as f:
        biop_args = dict([x.strip().split('=') for x in f.readlines()])

    return biop_args

def build_bioprospector_cmds(args):
    '''
    Build a list of commands to pass to subprocess to run bioprospector
    '''
    print("Building BioProspector commands...")

    cmds_list = [] # store all n cmds

    biop_args = parse_bioprospector_config(args.biop_config)

    # make output directory for this bioprospector run
    base = os.path.basename(args.seq_file).split('.')[0]
    biop_str = f"W{biop_args['W']}_w{biop_args['w']}_G{biop_args['G']}_g{biop_args['g']}_d{biop_args['d']}_a{biop_args['a']}_n{biop_args['n']}"
    ts = time.strftime("%s",time.gmtime())
    raw_dir = f"{base}_{biop_str}_BIOP_RAW_{ts}"
    raw_dir_path = os.path.join(args.outdir,raw_dir)
    
    mkdir_cmd = ['mkdir',raw_dir_path]
    
    # execute mkdir command and catch error
    p = Popen(mkdir_cmd, stdout=PIPE, stderr=PIPE)
    output, error = p.communicate()
    if p.returncode != 0: 
        raise ValueError(f'mkdir command failed:{error.decode("utf-8")}')
       
#    result = subprocess.run(cmds, stdout=subprocess.PIPE,stderr=subprocess.PIPE)

    # make a command string for each i in num_runs
    for i in range(args.num_runs):
        # args from config file
        cmds = [biop_args['biop_exe']] # executable
        cmds += ['-W',biop_args['W']]  # width of first block
        cmds += ['-w',biop_args['w']]  # width of second block
        cmds += ['-G',biop_args['G']]  # max spacer distance
        cmds += ['-g',biop_args['g']]  # min spacer distance
        cmds += ['-d',biop_args['d']]  # strand directions to check
        cmds += ['-a',biop_args['a']]  # do all inputs have a site
        cmds += ['-n',biop_args['n']]  # number of BioProspector internal iterations
        
        # input 
        cmds += ['-i',args.seq_file]   # input sequences
        
        # output file
        outf = os.path.join(raw_dir_path,f'biop_run{i}.txt')
        cmds += ['-o',outf]

        cmds_list.append(cmds)

    return cmds_list, raw_dir_path

def run_bioprospector(cmds_list):
    '''
    Given a list of command line args for bioprospector,
    execute the commands
    '''
    print("Executing BioProspector...")
    # get list of processes to execute each BioP command
    proc_list = [Popen(cmd, stdout=PIPE, stderr=PIPE) for cmd in cmds_list]
    # I think this for loop shoudl exectue the processes in parallel
    for i,proc in enumerate(proc_list):
        print(f"Run {i+1} of {len(proc_list)}")

        # execute bioprospector command and catch error
        #p = Popen(cmds, stdout=PIPE, stderr=PIPE)
        output, error = proc.communicate()
        if proc.returncode != 0: 
            raise ValueError(f'BioProspector execution failed:{error.decode("utf-8")}')


def tag_max(df):
    '''
    Given a df of a BioProspector summary file, for each gene, determine the 
    max block count. Return the df with the max blocks tagged
    '''
    sub_dfs = []
    for seq_name,sub_df in df.groupby('seq_name'):
        max_count = sub_df['block_count'].max()
        sub_df['is_max'] = sub_df.apply(
            lambda row: True if row['block_count'] == max_count else False,axis=1
        )
        sub_dfs.append(sub_df)
    
    return pd.concat(sub_dfs)

def summarize_biop_replicates(biop_summ_files):
    '''
    Given a list of `i` summary files created after parsing a group of `n`
    BioPropsector runs
    '''
    # for each summary file, load it into a df 
    dfs = [] # collect summary dfs
    for i,f in enumerate(biop_summ_files):
        df = pd.read_csv(f,sep='\t')
        df['settings'] = f'rep{i}'  # add replicate id
        dfs.append(tag_max(df))     # tag the max motif block for each gene

    # combine all replicate dfs into 1
    comb_df = pd.concat(dfs)

    selected_seqs = []

    # after combining replicates, find the max blocks for each seq name
    # and count the number of times a given motif was called 'max'
    for loc, sub_df in comb_df.groupby('seq_name'):
        # filter to only "is max" rows
        max_sub_df = sub_df[sub_df['is_max']]
        # count the number of iterations where this sequence was a "max row"
        vc = max_sub_df['seq_block'].value_counts(ascending=False).rename_axis('sequence').reset_index(name='count')
        # what was the max number of times a sequence was a max row
        max_agreement = vc['count'].values[0]
        # choose the sequence(s) that were most often chose as the max
        selected = vc[vc['count']==max_agreement]['sequence'].values

        selected_seqs.append((loc,selected))


    return selected_seqs
        
def write_final_selected_seqs(selected_seqs, outf):
    '''
    After summarizing BioP replicates, dump final selections to a file. Include ties
    '''
    with open(outf, 'w') as f:
        for (loc, seq_list) in selected_seqs:
            f.write(f">{loc}\n")
            # might have been ties... currently writes both seqs
            for seq in seq_list:
                f.write(f"{seq}\n")



# +-------------+
# | MAIN SCRIPT |
# +-------------+

def main():
    # +------+
    # | ARGS |
    # +------+
    parser = argparse.ArgumentParser(description='Run set of upstream regions through BioProspector to search for promoter motifs')
    # Required args
    parser.add_argument('seq_file', help='fasta file containing promoter regions to search for motifs')
    parser.add_argument('biop_config', help='path to BioProspector config file')
    parser.add_argument('outdir', help='Output directory where results are written')
    # Optional args
    parser.add_argument('-n', '--num_runs', type=int, default=200,help='Number of times to run BioProspector (default 50)')
    parser.add_argument('-k', '--top_k', type=int, default=3,help='Number of top motif picks to report for each locus')
    
    args = parser.parse_args()
    
    # build list of command to execute
    cmds_list, biop_raw_dir = build_bioprospector_cmds(args)

    # run BioProspector via subprocess
    run_bioprospector(cmds_list)

    # parse all the raw BioProspector output files, summarize 
    # motifs found, and save output summary and selection files
    selection_outf = f"{biop_raw_dir}_SELECTION.fa"
    summary_outf = f"{biop_raw_dir}_SUMMARY.tsv"
    mov_outf = f"{biop_raw_dir}_TOP_{args.top_k}_MOV.tsv"

    bu.compile_and_select(
        args.top_k,
        biop_raw_dir, 
        args.seq_file, 
        selection_outf,
        summary_outf,
        mov_outf)


    # compile results from all `i` replicates
    #final_selected_seqs = summarize_biop_replicates(biop_summ_files)
    
    # write final output
    # final_selection_outf = os.path.join(args.outdir,"FINAL_BIOP_SELECTIONS.fa")
    # write_final_selected_seqs(final_selected_seqs, final_selection_outf)

    print("Done!")
    

if __name__ == '__main__':
    main()