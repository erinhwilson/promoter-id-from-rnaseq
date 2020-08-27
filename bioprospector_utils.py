# bioprospector_utils.py

import logomaker
from os import listdir
from os.path import isfile, join
import pandas as pd
import numpy as np
import re



class SeqMotifMatch_2B:
    '''
    Class to hold info about a particular region of a sequence that was called
    a "match" to a given BioProspector motif
    '''
    def __init__(self,
                 motif_id,
                 seq_name, 
                 seq_len, 
                 site_num, 
                 block1pos, 
                 block2pos,
                 block1seq,
                 block2seq):
        
        self.motif_id = motif_id
        self.name = seq_name
        self.seq_len = int(seq_len)
        self.site_num = site_num
        self.block1seq = block1seq
        self.block2seq = block2seq

        if not block1pos.startswith("f"):
            raise ValueError(f"Parsing motifs not on the forward (f) not yet implemented. See {block1pos})")
        self.block1pos = int(block1pos.strip("f "))
        
        if not block2pos.startswith("f"):
            raise ValueError(f"Parsing motifs not on the forward (f) not yet implemented. See {block2pos})")
        self.block2pos = int(block2pos.strip("f "))
        
        # calculate spacer
        block1pos_end = self.block1pos + len(self.block1seq)
        self.spacer = self.block2pos - block1pos_end
        
        # calculate distance from end of seq
        self.block1_dist2end = 1+ self.seq_len - block1pos_end
        self.block2_dist2end = 1+ self.seq_len - (self.block2pos + len(self.block2seq))
        
        
    def pprint(self):
        print(f"Seq: {self.name}")
        print(f"Motif {self.motif_id} match instance #{self.site_num}")
        print(f"Length: {self.seq_len}, Block 1: {self.block1pos}, Block 2: {self.block2pos}")
        print(f"{self.block1seq} -- ({self.spacer}) -- {self.block2seq}")
        print(f"[ {self.block1_dist2end} ] -- ---- -- [ {self.block2_dist2end} ]")
        
    def quick_summ(self):
        return f"[{self.block1pos}]{self.block1seq} -- ({self.spacer}) -- [{self.block2pos}]{self.block2seq}[end-->{self.block2_dist2end}]"
        
    def extract_seq_block(self,seq_lookup):
        '''
        Given a sequence look up dictionary, find the sequence who's key
        matches the name of this object and return the whole slice of the 
        sequence between the first and second block
        '''
        
        # make sure the name of this sequence match is indeed in the lookup
        # dict provided from the bioprospector object
        assert(self.name in seq_lookup)
        
        full_seq = seq_lookup[self.name]
        left_coord = self.block1pos - 1
        right_coord = self.block2pos + len(self.block2seq) - 1
        
        seq_block = full_seq[left_coord:right_coord]
        return seq_block
        

class Motif_2B:
    '''
    Individual motif identified by BioProspector
    '''
    def __init__(self,motif_block):
        # parse the first line of the block into the 
        # ID and seq pattern
        first_line = motif_block[0]
        
        id_pattern = "Motif #(.*):"
        id_re = re.search(id_pattern,first_line)
        if not id_re:
            raise ValueError(f"No motif id found in {first_line}")
        self.motif_id = id_re.group(1)
        
        motif_pattern = "\((.*)/(.*), (.*)/(.*)\)"
        motif_re = re.search(motif_pattern, first_line)
        if not motif_re:
            raise ValueError(f"Motif pattern not found in {first_line}")
        self.block1 = motif_re.group(1)
        self.block1rev = motif_re.group(2)
        self.block2 = motif_re.group(3)
        self.block2rev = motif_re.group(4)
        
        # get the motif score and sites
        motif_block = motif_block[2:] # scoot down the motif block to the score line
        assert motif_block[0].startswith("Width")
        score_line = motif_block[0]
        score_pattern = "MotifScore (.*); Sites (.*)"
        score_re = re.search(score_pattern, score_line)
        if not score_re:
            raise ValueError(f"Score pattern not found in {score_line}")
        self.motif_score = score_re.group(1)
        self.num_sites = score_re.group(2)
        
        # skip down to the fasta alignments
        cur_line = motif_block[0]
        while not cur_line.startswith(">"):
            motif_block = motif_block[1:]
            cur_line = motif_block[0]
        
        seq_matches = []
        while not cur_line.startswith("*****"):
            # get gene info
            assert(cur_line.startswith(">"))
            header_pattern = ">(.*)\tlen (.*)\tsite #(.*)\t(.*)\t(.*)"
            header_re = re.search(header_pattern,cur_line)
            if not header_re:
                raise ValueError(f"Header pattern not found in {cur_line}")
                
            # get sequence match info
            seq_line = motif_block[1]
            seq_2B_pattern = "(.*) (.*)"
            seq_re = re.search(seq_2B_pattern,seq_line)
            if not seq_re:
                raise ValueError(f"Sequence match pattern not found in {seq_line}")
                
            # build a seq match object
            # motif_id, 
            # h1 seq_name, 
            # h2 seq_len, 
            # h3 site_num, 
            # h4 block1pos, 
            # h5 block2pos, 
            # s1 block1seq, 
            # s2 block2seq
            new_seq_match = SeqMotifMatch_2B(self.motif_id, 
                                             header_re.group(1),
                                             header_re.group(2),
                                             header_re.group(3),
                                             header_re.group(4),
                                             header_re.group(5),
                                             seq_re.group(1),
                                             seq_re.group(2)
                                            )
            seq_matches.append(new_seq_match)
            

            # scoot to next pair of lines line
            motif_block = motif_block[2:]
            cur_line = motif_block[0]
            
        # final list of all sequence matches for this motif
        self.seq_matches = seq_matches
        
        
    
    def pprint(self):
        print(f"Motif {self.motif_id}")
        print(f"Block 1: {self.block1} ({self.block1rev})")
        print(f"Block 2: {self.block2} ({self.block2rev})")
        print(f"Score: {self.motif_score}, Sites: {self.num_sites}")
        print()
        print(f"Number of Seq Matches: {len(self.seq_matches)}\n")
        for sm in self.seq_matches:
            sm.pprint()
            print()
            
            
    def build_pwm(self):
        '''
        Given all the sequences that match this motif, pile them into 
        a PWM
        '''
        oh = { # one-hot look up
            'A':0,
            'C':1,
            'G':2,
            'T':3
        }
        # initialize empty PWM
        block1_pwm = np.array([np.zeros(4) for x in range(len(self.seq_matches[0].block1seq))])
        block2_pwm = np.array([np.zeros(4) for x in range(len(self.seq_matches[0].block2seq))])

        # loop through all matches for this motif
        for sm in self.seq_matches:
            # block 1
            for i,base in enumerate(sm.block1seq):
                block1_pwm[i][oh[base]] += 1
                    
            # block 2
            for i,base in enumerate(sm.block2seq):
                block2_pwm[i][oh[base]] += 1             
        
        # convert to frequencies
        # block 1
        for i,pos in enumerate(block1_pwm):
            for j,base in enumerate(pos):
                # divide count at each pos by total seqs
                block1_pwm[i][j] = base/len(self.seq_matches)
        
        # block 2
        for i,pos in enumerate(block2_pwm):
            for j,base in enumerate(pos):
                # divide count at each pos by total seqs
                block2_pwm[i][j] = base/len(self.seq_matches)
        
        self.b1_pwm = block1_pwm
        self.b2_pwm = block2_pwm
        return block1_pwm, block2_pwm
    
    
    def view_motif(self):
        pwm1, pwm2 = self.build_pwm()
        df1 = pd.DataFrame(pwm1, columns=['A','C','G','T'])
        df2 = pd.DataFrame(pwm2, columns=['A','C','G','T'])

        # initialize 2 panel figure
        fig,(ax1,ax2) = plt.subplots(1,2,sharey=True,figsize=[10,2])
        
        # block 1
        logo1 = logomaker.transform_matrix(df1,from_type='probability',to_type='information')
        crp_logo1 = logomaker.Logo(logo1,ax=ax1)
        
        # block 2
        logo2 = logomaker.transform_matrix(df2,from_type='probability',to_type='information')
        crp_logo2 = logomaker.Logo(logo2,ax=ax2,)
        
        # labels
        df1_title = f"Block 1 ({len(self.seq_matches)} sites)"
        df2_title = f"Block 2 ({len(self.seq_matches)} sites)"
        ax1.set_title(df1_title)
        ax2.set_title(df2_title)
        
        fig.suptitle(f"Motif {self.motif_id} (score:{self.motif_score})", fontsize=20)
        fig.subplots_adjust(top=0.7) # shift main title up 
        plt.ylim(0,2)
        
        plt.show()

    def get_color(row):
        pal = sns.color_palette("hls", 4)
        if row['value'] == 1.0 and row['instance']=='1':
            return pal.as_hex()[0]
        elif row['value'] == 2.0 and row['instance']=='1':
            return pal.as_hex()[1]
        elif row['value'] == 1.0 and row['instance']=='2':
            return pal.as_hex()[2]
        elif row['value'] == 2.0 and row['instance']=='2':
            return pal.as_hex()[3]
        elif row['value'] == -1.0:
            return 'white'
        else:
            return 'black'    

    def visualize_motif_locations(self):
        '''
        TODO: This shoudl probably be a method of a Motif object?
        
        Given a motif found by BioProspector, use Altir to visualize
        the locations of this motif on the sequences that matched it.
        
        Returns a vertically concatted chart for all the sequence matches
        store in this motif's m.seq_matches variable.
        '''

        rows = []
        # for each promoter sequence that had a match to m
        for sm in self.seq_matches:
            row = []
            #print("Promoter Match", sm.name)
            row.append(sm.name)
            row.append(sm.seq_len)
            row.append(sm.site_num)
            row.append(sm.block1pos)
            row.append(sm.block2pos)
            row.append(sm.block1seq)
            row.append(sm.block2seq)

            rows.append(row)

        # convert seq match data into a dataframe for altair
        biop_df = pd.DataFrame(rows,
                               columns=['seq_name',
                                        'seq_len',
                                        'instance',
                                        'b1pos',
                                        'b2pos',
                                        'b1seq',
                                        'b2seq'])


        # chart initialization
        charts = []
        ticks = [int(x) for x in list(np.arange(-100,0,10))]
        
        # group by seq name (be sure to capture multiple instances of a 
        # matches of the motif to a sequence
        biop_dfg =  biop_df.groupby('seq_name')
        for name in biop_dfg.groups:
            df = biop_dfg.get_group(name)

            arrs = [] # array of mostly 0s except with the motif overlaps
            instance_arrs = [] # array of which match instance this mask is for
            b1seq_arrs = [] # copy of b1seq (same value in each slot, format for altair)
            b2seq_arrs = [] # copy of b2seq (same value in each slot, format for altair)
            pos_arrs = [] # position label for each slot in array (from -100:0)
            
            # for each row for this seq_name (may have multiple match instances)
            for i, row in df.iterrows():
                arr = np.zeros(row.seq_len)
                b1_s = row['b1pos']
                b1_e = row['b1pos'] + len(row['b1seq'])
                b1_seq = row['b1seq']

                b2_s = row['b2pos']
                b2_e = row['b2pos'] + len(row['b2seq'])
                b2_seq = row['b2seq']

                # make a mask over the sequence where the match overlaps in the promoter seq
                arr[b1_s:b1_e+1] = 1
                arr[b2_s:b2_e+1] = 2

                # provide padding if seq is shorter than 100
                if len(arr)<100:
                    # array of -1s
                    padding = np.zeros(100-len(arr)) + -1 
                    arr = np.concatenate((padding,arr))

                # put all the data into lists (later will make a df for altair)
                arrs.append(arr)
                instance_arrs.append([row['instance'] for x in arr])
                b1seq_arrs.append([b1_seq for x in arr])
                b2seq_arrs.append([b2_seq for x in arr])
                pos_arrs.append(np.arange(-100,0))
                


            # collect the array mask into a new df
            mask_df = pd.DataFrame(np.concatenate(arrs).ravel(), columns=['value'])
            mask_df['instance'] = np.concatenate(instance_arrs).ravel()
            mask_df['b1_seq'] = np.concatenate(b1seq_arrs).ravel()
            mask_df['b2_seq'] = np.concatenate(b2seq_arrs).ravel()
            mask_df['pos'] = np.concatenate(pos_arrs).ravel()
            mask_df['color'] = mask_df.apply(lambda row: get_color(row),axis=1)

            # build altair heatmap for a particular sequence
            chart = alt.Chart(
                mask_df,
                title=f"Motif loc in {name}"
            ).mark_rect().encode(
                x=alt.X('pos:O',axis=alt.Axis(values=ticks)),
                y=alt.Y('instance:O'),
                color=alt.Color('color:N',scale=None),
                tooltip=["b1_seq:N",'b2_seq:N'],
            ).properties(
                width=700)
            
            # collect charts in a list
            charts.append(chart)
            
        # vertically concat all charts for each motif
        mega_chart = charts[0]
        for c in charts[1:]:
            mega_chart = alt.vconcat(mega_chart, c)

        return mega_chart

    def view_motif_and_locs(self):
        self.view_motif()
        display(self.visualize_motif_locations())
        
    def extract_seq_match_blocks(self,seq_lookup):
        '''
        Given the set of sequences that match this 2-block motif,
        extract the fully sequence containing and between the motif
        match
        '''
        
        name2block = {}
        for sm in self.seq_matches:
            # start a dict for each new name we see    
            if sm.name not in name2block:
                name2block[sm.name] = {}
            # get the sequence block, site num, and actual SeqMatch object out
            block = sm.extract_seq_block(seq_lookup)
            name2block[sm.name][sm.site_num] = {
                'seq': block,
                'sm' : sm
            }
        
        return name2block


class BioProspectorResult:
    '''
    Given an output file from BioProspector, parse it into the 
    motifs and sequence match results contained within it
    '''
    def __init__(self,biop_filename, seq_filename):
        self.biop_filename = biop_filename
        self.seq_filename = seq_filename
        self.motifs = []
        
        motif_blocks = []
        with open(biop_filename,'r') as f:
            lines = [x.strip() for x in f.readlines()]
            cur = lines[0]
            # start with the first line and skip down past the "try" lines to the motif info
            while not cur.startswith("The highest scoring"):
                lines = lines[1:]
                cur = lines[0]

            # made it to the start of the motif info
            assert(lines[1].startswith("Motif #1"))
            # scoot down one more to the first motif
            lines = lines[1:]

            # break out the file into its raw sections of motifs.
            # specifically, collect the row index of lines that start
            # with "Motif #"
            motif_start_line = []
            for i,line in enumerate(lines):
                # stop when we reach total time:
                if line.startswith("Total time"):
                    break

                if line.startswith("Motif #"):
                    motif_start_line.append(i)

            # now I know which lines the new motifs start on. slice them out into individual blocks
            for i in range(len(motif_start_line)-1):
                motif_blocks.append(lines[motif_start_line[i]:motif_start_line[i+1]])
            # add the last block through the end
            motif_blocks.append(lines[motif_start_line[-1]:]) 

        
        # with the motif blocks of raw Bioprospector output collected,
        # turn each block into a 2-Block motif object
        for mb in motif_blocks:
            self.motifs.append(Motif_2B(mb))
            
        # load the sequence (promoter) file and store a lookup dict between the 
        # full promoter name and the sequence
        promoters = load_promoter_seqs(self.seq_filename)
        self.name2seq_lookup = dict([(x[1],x[2]) for x in promoters])
        
        #print(f"Done parsing BioProspector file {self.biop_filename}")
        #print(f"Found {len(self.motifs)} motifs")
        
    def view_motifs(self):
        for m in self.motifs:
            m.view_motif()
            
    def view_motifs_and_locs(self):
        for m in self.motifs:
            m.view_motif_and_locs()
    
    def extract_sequence_blocks(self):
        '''
        For all the motifs found this bioprospector run, extract all the 
        the seq blocks
        '''
        motif_seq2block = {}
        for m in self.motifs:
            seq2block = m.extract_seq_match_blocks(self.name2seq_lookup)
            motif_seq2block[m.motif_id] = seq2block
        
        return motif_seq2block

def load_promoter_seqs(filename):
    '''
    Load fasta file of promoters into ID, desc, and seq
    '''
    proms = []
    with open(filename,'r') as f:
        for line in f:
            if line.startswith(">"):
                full_header = line.strip()[1:]
                locus_tag = full_header.split('|')[0]
            else:
                seq = line.strip()
                proms.append((locus_tag,full_header,seq))
                
    return proms


def compile_all_biop_results(biop_files, promoter_file,top_motif_only=False, verbose=False):
    '''
    Given a list of BioProspector output files and a list of 
    promoters use as input to those BioProspector files (promoters
    must be the same inputs to each BioP file! Otherwise the BioP
    file is reporting matches to different sequences/sections of
    sequences)
    '''
    promoters = load_promoter_seqs(promoter_file)
    
    rows = []
    # loop through promoters
    for lt,name,seq in promoters:
        row = []
        if verbose:
            print(f"\n> {name}")
        # for each bioP output file
        for biop_file in biop_files:
            # load the BioP result
            biop_result = BioProspectorResult(biop_file,promoter_file)
            # extract the full motif sequence blocks from the BioP match coordinates
            bpr_blocks = biop_result.extract_sequence_blocks()
            # for each motif in this BioP result
            for motif_id in bpr_blocks:
                # if in top_motif_only mode, skip the remaining motifs
                if top_motif_only and int(motif_id) > 1:
                    continue
                else:
                    # skip if this promoter didn't have a match to this motif
                    if name not in bpr_blocks[motif_id]:
                        if verbose:
                            print(f"{name} not in {biop_file}, motif {motif_id}")
                    # otherwise, extract the motif info
                    else:
                        # loop through multiple instances
                        for instance in bpr_blocks[motif_id][name]:
                            
                            sm_seq = bpr_blocks[motif_id][name][instance]['seq']
                            sm_obj = bpr_blocks[motif_id][name][instance]['sm']
                            if verbose:
                                print(f"{sm_seq} {sm_obj.quick_summ()} (biop {biop_file}, motif {motif_id}, inst {instance})")

                            row = [name,sm_seq,sm_obj.quick_summ(),biop_file,motif_id,instance]
                            rows.append(row)
    
    df = pd.DataFrame(rows, columns=['seq_name','seq_block','block_summ','biop_file','motif_id','instance'])
    return df



def select_motifs_with_most_agreement(df,
                                      seq_selection_out,
                                      summary_out,
                                      verbose=False):
    '''
    Given a df compiling the motif match results from several bioP runs,
    figure out which motif for each sequence had the most agreement between
    runs of the BioP program.
    
    Output a fasta file of the selected sequences and a summary file for reference.
    '''
    selection_rows = []
    summary_rows = []
    
    # group the compiled bioP df by promoter sequence
    dfg = df.groupby('seq_name')

    # for each sequence, count the number of times each seq match was 
    # identified by BioP
    for seq_name,sub_df in dfg:
        #sub_df = dfg.get_group(seq_name)
        if verbose:
            print(f"Getting seqs for {seq_name}")
        # use value counts function to count number of occurences of each 
        # block summary. Then some nasty reformatting/column renaming
        vc = sub_df['block_summ'].value_counts(ascending=False)
        vc = vc.reset_index().rename(columns={'block_summ':'count','index':'block_summ'})

        # Now we know the order of seq blocks from high to low. 
        # Compile a summary file and final seq selection file
        
        # group our sub_df by block summaries
        dfgg = sub_df.groupby('block_summ')

        # now that we have an ordering of the most frequent matches, print out a summary
        for i,row in vc.iterrows():
            
            block_summ = row['block_summ']
            count = row['count']
            # use this block summ to get the group from the sub_df
            block_df = dfgg.get_group(block_summ)
            assert(len(block_df['seq_block'].unique())==1)
            seq_block = block_df['seq_block'].values[0]
            agreements = [f"{bio_desc}_MOTIF-{motif_id}-{inst}" for (bio_desc,motif_id,inst) in block_df[['biop_file','motif_id','instance']].values]
            
            # if this the block with the most agreement, add it to the selection fasta
            if i == 0:
                if verbose:
                    print(f"Top selection:{block_summ} ({count})")
                selection_rows.append(f"> {seq_name}\n")
                selection_rows.append(f"{seq_block}\n")
                
            # regardless of whether its the best match, write out to the summary file
            summary_rows.append(f"{seq_name}\t{count}\t{seq_block}\t{block_summ}\t{agreements}\n")
    
    # write out the selection and summary
    with open(seq_selection_out,'w') as f:
        for row in selection_rows:
            f.write(row)

    with open(summary_out,'w') as f:
        f.write(f"seq_name\tblock_count\tseq_block\tblock_summ\tagreements\n")
        for row in summary_rows:
            f.write(row)
            
    print("\nDone writing selection and summary to:")
    print("Selection:",seq_selection_out)
    print("Summary:",summary_out)

def write_margin_of_victory_results(top_k, summary_outf, margin_of_victory_outf):
    '''
    After having written out the bioprospector summary file, 
    use that to determine the top k motifs with the most
    "BioProspector Votes" and calculate the margin of victory 
    for each motif, for each locus
    '''
    # load up the summary file 
    sum_df = pd.read_csv(summary_outf,sep='\t')

    rows = []
    for loc, df in sum_df.groupby('seq_name'):
        # get the top k voted motifs
        top_df = df.sort_values('block_count', ascending=False).head(top_k+1)
        
        for i in range(top_k):
            seq, block_summ, block_count = top_df[['seq_block', 'block_summ', 'block_count']].values[i]
            next_block_count = top_df['block_count'].values[i+1]
            # margin of victory (difference between current and next most popular motif)
            mov = block_count - next_block_count

            # add row to data
            row = [loc,seq,mov,block_summ, block_count]
            rows.append(row)
            
    mov_df = pd.DataFrame(rows, columns=['loc','sequence','margin_of_victory','motif_summ','raw_votes'])
    mov_df.to_csv(margin_of_victory_outf, sep='\t',index=False)

    print("Margin of victory:",margin_of_victory_outf)

def compile_and_select(top_k,
                       biop_raw_dir, 
                       promoter_file, 
                       selection_outf, 
                       summary_outf, 
                       mov_outf,
                       top_motif_only=False):
    
    # get list of files from bioprospector output directory
    biop_files = [join(biop_raw_dir,f) for f in listdir(biop_raw_dir) if isfile(join(biop_raw_dir, f))]

    # compile results
    print("Compiling BioProspector result files...")
    df = compile_all_biop_results(biop_files, 
                                  promoter_file,
                                  top_motif_only=top_motif_only)
    
    # select the top motif for each seq
    print("Selecting best motif for each sequence...")
    select_motifs_with_most_agreement(df,
                                      selection_outf,
                                      summary_outf,
                                      verbose=True)

    # summarize the margin of error
    print("Reporting top motifs with margins of victory...")
    write_margin_of_victory_results(top_k, summary_outf, mov_outf)

    

