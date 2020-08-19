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
from subprocess import Popen, PIPE
import time


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
    cmds_list = [] # store all n cmds

    biop_args = parse_bioprospector_config(args.biop_config)
    print(biop_args)

    # make output directory for this bioprospector run
    base = os.path.basename(args.seq_file).split('.')[0]
    biop_str = f"W{biop_args['W']}_w{biop_args['w']}_G{biop_args['G']}_g{biop_args['g']}_d{biop_args['d']}_a{biop_args['a']}_n{biop_args['n']}"
    ts = time.strftime("%s",time.gmtime())
    raw_dir = f"BIOP_RAW_{base}_{biop_str}_{ts}"
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
    print("Executing BioProspector")
    for i,cmds in enumerate(cmds_list):
        print(f"Run {i+1} of {len(cmds_list)}")

        # execute bioprospector command and catch error
        p = Popen(cmds, stdout=PIPE, stderr=PIPE)
        output, error = p.communicate()
        if p.returncode != 0: 
            raise ValueError(f'BioProspector execution failed:{error.decode("utf-8")}')
      



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
    parser.add_argument('-n', '--num_runs', type=int, default=2,help='Number of times to run BioProspector (default 20)')
    
    args = parser.parse_args()
    print(args.seq_file)
    print(args.biop_config)
    print(args.outdir)
    print(args.num_runs)

    cmds_list, biop_raw_dir = build_bioprospector_cmds(args)
    
    run_bioprospector(cmds_list)

    # +------+
    # | SAVE |
    # +------+
    # save to a fasta file
    # out_path = write_fasta_file(args, loci, upstream_regions)    

    # print(f"Output written to {out_path}")
    # print("Done!")
    

if __name__ == '__main__':
    main()