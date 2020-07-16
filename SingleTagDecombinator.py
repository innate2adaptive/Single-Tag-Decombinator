##################
##### UPDATE #####
##################

# This script has been updated to work for very short read single cell data. It will search for TCRs by looking for 
# a single tag within a read, unlike Classic Decombinator which searches for both a V and J tag per read. 
# Thomas Peacock, February 2017, UCL

# Classic Decombinator written by 
# James M. Heather, August 2016, UCL
# https://innate2adaptive.github.io/Decombinator/


##################
### BACKGROUND ###
##################

# Searches FASTQ reads (produced through Demultiplexor.py) for rearranged TCR chains
# Can currently analyse human and mouse TCRs, both alpha/beta and gamma/delta chains
  # NB: Human alpha/beta TCRs are the most thoroughly tested, due to the nature of the data we generated. YMMV.

# Current version (v3) is optimised for interpretation of data generated using our wet lab protocol, but could be modified to work on any data.

# Script represents improvements upon a combination of the two previously in use Decombinator versions
  # i.e. Decombinator V2.2 (written primarily by Nic Thomas, see Thomas et al, Bioinformatics 2013, DOI: 10.1093/bioinformatics/btt004)
  # and vDCR (which was v1.4 modified by James Heather, see Heather et al, Frontiers in Immunology 2016, DOI: 10.3389/fimmu.2015.00644)
  # Now faster, more accurate and easier to use than either of the previous versions.
  
##################
###### INPUT #####
##################

# As with entire pipeline, Decombintator is run using command line arguments to provide user parameters
  # All arguments can be read by viewing the help data, by running python Decombintator.py -h

# Takes FASTQ reads produced by Demultiplexor.py (unzipped or gzipped), which is the minimum required command line input, using the -fq flag
  # NB: Data must have been generated using the appropriate 5'RACE ligation protocol, using the correct SP2-I8-6N-I8-6N oligonucleotide

# The TCR chain locus to look for can be explicitly specified using the -c flag 
  # Users can use their choice of chain identifiers from this list (case insensitive): a/b/g/d/alpha/beta/gamma/delta/TRA/TRB/TRG/TRD/TCRA/TCRB/TCRG/TCRD
  # If no chain is provided (or if users which to minimise input arguments), script can infer chain from the FASTQ filename
    # I.e. "alpha_sample.fq" would be searched for alpha chain recombinations
    # NB: This autodetection only works if there is only ONE TCR locus present in name (which must be spelt out in full)

# Other optional flags:
  
  # -s/--supresssummary: Supress the production of a summary file containing details of the run into a 'Logs' directory. 
      
  # -dz/--dontgzip: Suppress the automatic compression of output demultiplexed FASTQ files with gzip. 
    # Using this flag makes the script execute faster, but data will require more storage space. 
    
  # -dc/--dontcount: Suppress the whether or not to show the running line count, every 100,000 reads. 
    # Helps in monitoring progress of large batches.
  
  # -dk/--dontcheck: Suppress the FASTQ sanity check. 
    # Strongly recommended to leave alone: sanity check inspects first FASTQ read for basic FASTQ parameters.
  
  # -pf/--prefix: Allows users to specify the prefix of the Decombinator TCR index files produced. Default = 'dcr_'
  
  # -ex/--extension: Allows users to specify the file extension of the Decombinator TCR index files produced. Default = '.n12'

  # -or/--orientation: Allows users to specify which DNA orientations to check for TCR reads. Default = reverse only, as that's what the protocol produces.
    # This will likely need to be changed for analysing data produced by protocols other than our own.

  # -tg/--tags: Allows users to specify which tag set they wish to use. For human alpha/beta TCRs, a new 'extended' tag set is recommended, as it covers more genes.
    # Unfortunately an extended tag set is only currently available for human a/b genes.

  # -sp/--species: Current options are only human or mouse. Help could potentially be provided for generation of tags for different species upon request.
  
  # -N/--allowNs: Provides users the option to allow 'N's (ambiguous base calls), overriding the filter that typically removes rearrangements that contain them.
    # Users are recommended to not allow Ns, as such bases are both themselves low quality data and predict reads that are generally less trustworthy.
    
  # -ln/--lenthreshold: The length threshold which (the inter-tag region of) successful rearrangements must be under to be accepted. Default = 130.
  
  # -tfdir/--tagfastadir: The path to a local copy of a folder containing the FASTA and Decombinator tag files required for offline analysis.
    # Ordinarily such files can be downloaded on the fly, reducing local clutter.
    # By default the script looks for the required files in the present working directory, then in a subdirectory called "Decombinator-Tags-FASTAs", then online.
    # Files are hosted on GitHub, here: https://github.com/innate2adaptive/Decombinator-Tags-FASTAs

  # -nbc/--nobarcoding: Run Decombinator without any barcoding, i.e. use the whole read. 
    # Recommended when running on data not produced using the Innate2Adaptive lab's ligation-mediated amplification protocol

##################
##### OUTPUT #####  
##################

# Produces a '.n12' file by default, which is a standard comma-delimited Decombinator output file with several additional fields:
  # V index, J index, # V deletions, # J deletions, insert, ID, TCR sequence, TCR quality, barcode sequence, barcode quality
  # NB The TCR sequence given here is the 'inter-tag' region, i.e. the sequence between the start of the found V tag the end of the found J tag 

##################
#### PACKAGES ####  
##################

from __future__ import division
import sys          
import os
import urllib2
import string
import collections as coll
import argparse
import gzip
import Levenshtein as lev
import collections
from Bio import SeqIO
from Bio.Seq import Seq
from acora import AcoraBuilder
from time import time, strftime

__version__ = '3.1'

##########################################################
############# READ IN COMMAND LINE ARGUMENTS #############
##########################################################

def args():
  """args(): Obtains command line arguments which dictate the script's behaviour"""

  # Help flag
  parser = argparse.ArgumentParser(
      description='Decombinator v3.1: find rearranged TCR sequences in HTS data. Please go to https://innate2adaptive.github.io/Decombinator/ for more details.')
  # Add arguments
  parser.add_argument(
      '-fq', '--fastq', type=str, help='Correctly demultiplexed/processed FASTQ file containing TCR reads', required=True)
  parser.add_argument(
      '-fq2', '--fastq2', type=str, help='Optional second file to be to be merged with original file before analysis (e.g. for merging multiple'\
      +'reads R1 and R2 for simultaneous TCR search)', required=False)
  parser.add_argument(
      '-c', '--chain', type=str, help='TCR chain (a/b/g/d)', required=False)
  parser.add_argument(
      '-s', '--suppresssummary', action='store_true', help='Suppress the production of summary data log file', required=False)
  parser.add_argument(
      '-dz', '--dontgzip', action='store_true', help='Stop the output FASTQ files automatically being compressed with gzip', required=False)
  parser.add_argument(
      '-dk', '--dontcheck', action='store_true', help='Skip the FASTQ check', required=False, default=False)  
  parser.add_argument(
      '-dc', '--dontcount', action='store_true', help='Stop Decombinator printing a running count', required=False)
  parser.add_argument(
      '-ex', '--extension', type=str, help='Specify the file extension of the output DCR file. Default = \"n12\"', required=False, default="n12")
  parser.add_argument(
      '-pf', '--prefix', type=str, help='Specify the prefix of the output DCR file. Default = \"dcr_\"', required=False, default="dcr_")
  parser.add_argument(
      '-or', '--orientation', type=str, help='Specify the orientation to search in (forward/reverse/either/both). Default = reverse', required=False, default="reverse")  
  parser.add_argument(
      '-tg', '--tags', type=str, help='Specify which Decombinator tag set to use (extended or original). Default = extended', required=False, default="extended")
  parser.add_argument(
      '-sp', '--species', type=str, help='Specify which species TCR repertoire the data consists of (human or mouse). Default = human', required=False, default="human")
  parser.add_argument(
      '-N', '--allowNs', action='store_true', help='Whether to allow VJ rearrangements containing ambiguous base calls (\'N\'). Default = False', required=False)
  parser.add_argument(
      '-ln', '--lenthreshold', type=int, help='Acceptable threshold for inter-tag (V to J) sequence length. Default = 130', required=False, default=130)
  parser.add_argument(
      '-tfdir', '--tagfastadir', type=str, help='Path to folder containing TCR FASTA and Decombinator tag files, for offline analysis.', \
      required=False, default="Decombinator-Tags-FASTAs")
  parser.add_argument(
      '-nbc', '--nobarcoding', action='store_true', help='Option to run Decombinator without barcoding, i.e. so as to run on data produced by any protocol.', required=False)
  parser.add_argument(
      '-np', '--nproc', type=int, help='Number of cores for multprocessing alignment', required=False, default=None)
  return parser.parse_args()

##########################################################
############# FASTQ SANITY CHECK AND PARSING #############
##########################################################

def fastq_check(infile):
  """fastq_check(file): Performs a rudimentary sanity check to see whether a file is indeed a FASTQ file"""
  
  success = True
    
  #if infile.endswith('.gz'):
  with opener(infile) as possfq:
    try:
      read = [next(possfq) for x in range(4)]
    except:
      print "There are fewer than four lines in this file, and thus it is not a valid FASTQ file. Please check input and try again."
      sys.exit()
  # @ check
  if read[0][0] <> "@":
    success = False
  # Descriptor check
  if read[2][0] <> "+":
    success = False
  # Read/quality match check
  if len(read[1]) <> len(read[3]):
    success = False
  return(success)

def revcomp(read):
  """rc(read): Wrapper for SeqIO reverse complement function"""
  return str(Seq(read).reverse_complement())

def read_tcr_file(species, tagset, chain, gene, filetype, expected_dir_name):
  """ Reads in the FASTA and tag data for the appropriate TCR locus """
  
  # Define expected file name
  expected_file = species + "_" + tagset + "_" + "TR" + chain.upper() + gene.upper() + "." + filetype
  # expected_file = species + "_" + tagset + "_" + "TR" + chain.upper() + gene.upper() + "." + filetype

  # First check whether the files are available locally (in pwd or in bundled directory)
  if os.path.isfile(expected_file):
    fl = expected_file
    fl_opener = open
  elif os.path.isfile(expected_dir_name + os.sep + expected_file):
    fl = expected_dir_name + os.sep + expected_file
    fl_opener = open
  else:
    try:
      fl = "https://raw.githubusercontent.com/innate2adaptive/Decombinator-Tags-FASTAs/master/" + expected_file
      urllib2.urlopen(urllib2.Request(fl))      # Request URL, see whether is found
      fl_opener = urllib2.urlopen
    except:
      print "Cannot find following file locally or online:", expected_file
      print "Please either run Decombinator with internet access, or point Decombinator to local copies of the tag and FASTA files with the \'-tf\' flag."
      sys.exit()
  
  # Return opened file, for either FASTA or tag file parsing
  return fl_opener(fl)

def readfq(fp): 
    """
    readfq(file):Heng Li's Python implementation of his readfq function 
    https://github.com/lh3/readfq/blob/master/readfq.py
    """
    
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

#####################################
############# DECOMBINE #############
#####################################

def vanalysis(read):

  hold_v = v_key.findall(read)
  
  if hold_v:
    if len(hold_v) > 1:
      counts['multiple_v_matches'] += 1
      return
    v_match = v_seqs.index(hold_v[0][0]) # Assigns V
    temp_end_v = hold_v[0][1] + jump_to_end_v[v_match] - 1 # Finds where the end of a full V would be
    
    v_seq_start = hold_v[0][1]      
    end_v_v_dels = get_v_deletions( read, v_match, temp_end_v, v_regions )      
    if end_v_v_dels: # If the number of deletions has been found
      return v_match, end_v_v_dels[0], end_v_v_dels[1], v_seq_start
      
  else:
    
    hold_v1 = half1_v_key.findall(read)
    
    if hold_v1:
      for i in range(len(hold_v1)):
        indices = [y for y, x in enumerate(half1_v_seqs) if x == hold_v1[i][0] ]
        for k in indices:
          if len(v_seqs[k]) == len(read[hold_v1[i][1]:hold_v1[i][1]+len(v_seqs[half1_v_seqs.index(hold_v1[i][0])])]):
            if lev.hamming( v_seqs[k], read[hold_v1[i][1]:hold_v1[i][1]+len(v_seqs[k])] ) <= 1:
              counts['verr2'] += 1
              v_match = k
              temp_end_v = hold_v1[i][1] + jump_to_end_v[v_match] - 1 # Finds where the end of a full V would be
              end_v_v_dels = get_v_deletions( read, v_match, temp_end_v, v_regions )
              if end_v_v_dels:
                v_seq_start = hold_v1[i][1]  
                return v_match, end_v_v_dels[0], end_v_v_dels[1], v_seq_start
      counts['foundv1notv2'] += 1
      return
    
    else:
      
      hold_v2 = half2_v_key.findall(read)
      if hold_v2:
        for i in range(len(hold_v2)):
          indices = [y for y, x in enumerate(half2_v_seqs) if x == hold_v2[i][0] ]
          for k in indices:
            if len(v_seqs[k]) == len(read[hold_v2[i][1]-v_half_split:hold_v2[i][1]-v_half_split+len(v_seqs[half2_v_seqs.index(hold_v2[i][0])])]):
              if lev.hamming( v_seqs[k], read[hold_v2[i][1]-v_half_split:hold_v2[i][1]+len(v_seqs[k])-v_half_split] ) <= 1:
                counts['verr1'] += 1
                v_match = k
                temp_end_v = hold_v2[i][1] + jump_to_end_v[v_match] - v_half_split - 1 # Finds where the end of a full V would be
                end_v_v_dels = get_v_deletions( read, v_match, temp_end_v, v_regions )
                if end_v_v_dels:
                  v_seq_start = hold_v2[i][1] - v_half_split      
                  return v_match, end_v_v_dels[0], end_v_v_dels[1], v_seq_start
        counts['foundv2notv1'] += 1
        return
              
      else:
        counts['no_vtags_found'] += 1
        return
      
def janalysis(read):
  
  hold_j = j_key.findall(read)
  
  if hold_j:
    if len(hold_j) > 1:
      counts['multiple_j_matches'] += 1
      return
  
    j_match = j_seqs.index(hold_j[0][0]) # Assigns J
    temp_start_j = hold_j[0][1] - jump_to_start_j[j_match] # Finds where the start of a full J would be
    
    j_seq_end = hold_j[0][1] + len(hold_j[0][0])      
        
    start_j_j_dels = get_j_deletions( read, j_match, temp_start_j, j_regions )
    
    if start_j_j_dels: # If the number of deletions has been found

      return j_match, start_j_j_dels[0], start_j_j_dels[1], j_seq_end
          
  else:
    
    hold_j1 = half1_j_key.findall(read)
    if hold_j1:
      for i in range(len(hold_j1)):
        indices = [y for y, x in enumerate(half1_j_seqs) if x == hold_j1[i][0] ]
        for k in indices:
          if len(j_seqs[k]) == len(read[hold_j1[i][1]:hold_j1[i][1]+len(j_seqs[half1_j_seqs.index(hold_j1[i][0])])]):
            if lev.hamming( j_seqs[k], read[hold_j1[i][1]:hold_j1[i][1]+len(j_seqs[k])] ) <= 1:
              counts['jerr2'] += 1
              j_match = k
              temp_start_j = hold_j1[i][1] - jump_to_start_j[j_match] # Finds where the start of a full J would be
              j_seq_end = hold_j1[i][1] + len(hold_j1[i][0]) + j_half_split                                              
              start_j_j_dels = get_j_deletions( read, j_match, temp_start_j, j_regions )
              if start_j_j_dels:
                return j_match, start_j_j_dels[0], start_j_j_dels[1], j_seq_end
      counts['foundj1notj2'] += 1
      return              
            
    else:        
      hold_j2 = half2_j_key.findall(read)
      if hold_j2:
        for i in range(len(hold_j2)):
          indices = [y for y, x in enumerate(half2_j_seqs) if x == hold_j2[i][0] ]
          for k in indices:
            if len(j_seqs[k]) == len(read[hold_j2[i][1]-j_half_split:hold_j2[i][1]-j_half_split+len(j_seqs[half2_j_seqs.index(hold_j2[i][0])])]):
              if lev.hamming( j_seqs[k], read[hold_j2[i][1]-j_half_split:hold_j2[i][1]+len(j_seqs[k])-j_half_split] ) <= 1:
                counts['jerr1'] += 1
                j_match = k
                temp_start_j = hold_j2[i][1] - jump_to_start_j[j_match] - j_half_split # Finds where the start of a full J would be
                j_seq_end = hold_j2[i][1] + len(hold_j2[i][0])                                                
                start_j_j_dels = get_j_deletions( read, j_match, temp_start_j, j_regions )
                if start_j_j_dels:
                  return j_match, start_j_j_dels[0], start_j_j_dels[1], j_seq_end
        counts['foundv2notv1'] += 1
        return
      
      else:
         counts['no_j_assigned'] += 1
         return
       
def dcr(read, inputargs, chain_order):

  """dcr(read): Core function which checks a read (in the given frame) for a rearranged TCR of the specified chain.
    Returns a list giving: V gene index (if found), J gene index (if found), seq from end of V tag to end or read
    (or from start of read to start of J tag), position of end of V tag in read (or position of start of read), 
    position of end of read (or position of start of J tag in read). The last two fields are used to find the 
    appropriate quality score of the relevant sequence.
     """

  vdat = vanalysis(read)
  
  jdat = janalysis(read)

  if vdat:
    v_chain_order = []
    lim = 0
    v_seq_count = 0
    for i in range(len(chain_order)):
      if 'v' in chain_order[i]:
        v_chain_order.append(chain_order[i])

    for i in range(len(v_chain_order)):
      v_seq_count = lim
      lim += v_chain_order[i][2]
      if vdat[0] < lim:
        chain_type = v_chain_order[i][0]
        break
    vindex = vdat[0] - v_seq_count #adjusts to give correct index if analysing for multiple chainnams


  if jdat:
    j_chain_order = []
    lim = 0
    j_seq_count = 0
    for i in range(len(chain_order)):
      if 'j' in chain_order[i]:
        j_chain_order.append(chain_order[i])

    for i in range(len(j_chain_order)):
      j_seq_count = lim
      lim += j_chain_order[i][2]
      if jdat[0] < lim:
        chain_type = j_chain_order[i][0]
        break
    jindex = jdat[0] - j_seq_count

  if not vdat:
    #vdat = ["n/a"]
    vindex = "n/a"
  if not jdat:
    #jdat = ["n/a"]
    jindex = "n/a"

  if jindex != "n/a":
    start_of_j = jdat[3]-len(j_seqs[jdat[0]])   
    j_details = [vindex, jindex, read[0:jdat[3]], 0, jdat[3], chain_type]
    return j_details

  elif vindex != "n/a":
    end_of_v = vdat[3]+len(v_seqs[vdat[0]])
    v_details = [vindex, jindex, read[vdat[3]:len(read)], vdat[3], len(read), chain_type]
    return v_details
  else:
    counts['VJ_assignment_failed'] += 1
    return
  
###########################################################
############# ANCILLARY DECOMBINING FUNCTIONS #############
###########################################################

def get_chain(inputargs):

  nochain_error = "TCR chains not recognised. \n \
      Please either include at least one chain name in the file name (i.e. alpha/beta/gamma/delta),\n \
      or use the \'-c\' flag with an explicit chain option (a/b/g/d, case-insensitive).\n \
      Decombinator accepts a/b/g/d chains only."

  global chainnams, counts
  counts = coll.Counter()

  chainnams = {"a": "alpha", "b": "beta", "g": "gamma", "d": "delta"}
   
  # Detect whether chain specified in filename
  inner_filename_chains = [x for x in chainnams.values() if x in inputargs['fastq'].lower()]
  counts['chain_detected'] = len(inner_filename_chains)
  

  if inputargs['chain']:
    input_chains = inputargs['chain'].upper().split(" ")
  else:
    # If no chain provided, try and infer from filename
    if counts['chain_detected'] >= 1:
      
      input_chains = [c[0].upper() for c in inner_filename_chains]
    
    else:
      print nochain_error
      sys.exit()    
  chains = []

  if not input_chains:
    print nochain_error
    sys.exit()

  for i in range(len(input_chains)):
    if input_chains[i] in ['A', 'ALPHA', 'TRA', 'TCRA']:
      chains.append("a") 
    elif input_chains[i] in ['B', 'BETA', 'TRB', 'TCRB']:
      chains.append("b") 
    elif input_chains[i] in ['G', 'GAMMA', 'TRG', 'TCRG']:
      chains.append("g") 
    elif input_chains[i] in ['D', 'DELTA', 'TRD', 'TCRD']:
      chains.append("d") 
    else:
      print nochain_error
      sys.exit()
    chains = list(set(chains))
  return chains

def flatten(l):
  return [item for sublist in l for item in sublist]



def import_tcr_info(inputargs):
  """ import_tcr_info: Gathers the required TCR chain information for Decombining """
    
  # Get chain information
  global chain

  chain = get_chain(inputargs)

  #################################################
  ############# GET GENES, BUILD TRIE #############
  #################################################

  print 'Importing TCR', ", ".join(map(chainnams.__getitem__, chain)), 'gene sequences...'

  # First check that valid tag/species combinations have been used
  if inputargs['tags'] == "extended" and inputargs['species'] == "mouse":
    print "Please note that there is currently no extended tag set for mouse TCR genes.\n \
    Decombinator will now switch the tag set in use from \'extended\' to \'original\'.\n \
    In future, consider editing the script to change the default, or use the appropriate flags (-sp mouse -tg original)."
    inputargs['tags'] = "original"
  
  if inputargs['tags'] == "extended" and ( 'g' in chain or 'd' in chain ):
  
    print "Please note that there is currently no extended tag set for gamma/delta TCR genes.\n \
    Decombinator will now switch the tag set in use from \'extended\' to \'original\' for these chains.\n \
    In future, consider editing the script to change the default, or use the appropriate flags."
    inputargs['tags'] = "original"

  # Set tag split position, and check tag set. Note that original tags use shorter length J half tags, as these tags were originally shorter.

  global v_half_split, j_half_split
  if inputargs['tags'] == "extended":
    v_half_split, j_half_split = [10,10] 
  elif inputargs['tags'] == "original":
    v_half_split, j_half_split = [10,6] 
  else:
    print "Tag set unrecognised; should be either \'extended\' or \'original\' for human, or just \'original\' for mouse. \n \
    Please check tag set and species flag."
    sys.exit()

    
  # Check species information
  if inputargs['species'] not in ["human", "mouse"]:
    print "Species not recognised. Please select either \'human\' (default) or \'mouse\'.\n \
    If mouse is required by default, consider changing the default value in the script."
    sys.exit()    
    
  # Look for tag and V/J fasta and tag files: if these cannot be found in the working directory, source them from GitHub repositories
    # Note that fasta/tag files fit the pattern "species_tagset_gene.[fasta/tags]"
    # I.e. "[human/mouse]_[extended/original]_TR[A/B/G/D][V/J].[fasta/tags]"
  
  chain_order= []

  for gene in ['v', 'j']:

    # Get FASTA data
    fasta_holder = []

    for i in range(len(chain)):
      fasta_file = read_tcr_file(inputargs['species'], inputargs['tags'], chain[i], gene, "fasta", inputargs['tagfastadir'])  
      fasta_holder.append(list(SeqIO.parse(fasta_file, "fasta")))
      fasta_file.close()
      chain
    globals()[gene + "_genes"] = flatten(fasta_holder)
    
    
    
    globals()[gene+"_regions"] = []
    for g in range(0, len(globals()[gene+"_genes"])):
        globals()[gene+"_regions"].append(string.upper(globals()[gene+"_genes"][g].seq))  
        
    # Get tag data

    gene_seq_holder = []  #initialise arrays
    half1_gene_seq_holder = []
    half2_gene_seq_holder = []
    jumpfunction_holder = []

    for i in range(len(chain)):
      tag_file = read_tcr_file(inputargs['species'], inputargs['tags'], chain[i], gene, "tags", inputargs['tagfastadir'])  # get tag data
      if gene == 'v': jumpfunction = "jump_to_end_v"
      elif gene == 'j': jumpfunction = "jump_to_start_j"
      tag_info_holder = globals()["get_"+gene+"_tags"](tag_file, globals()[gene+"_half_split"])
      gene_seq_holder.append(tag_info_holder[0])
      half1_gene_seq_holder.append(tag_info_holder[1])
      half2_gene_seq_holder.append(tag_info_holder[2])
      jumpfunction_holder.append(tag_info_holder[3])
      chain_order.append([chain[i],gene, len(gene_seq_holder[i])])
      tag_file.close()

    globals()[gene+"_seqs"] = flatten(gene_seq_holder)
    globals()["half1_"+gene+"_seqs"] = flatten(half1_gene_seq_holder)
    globals()["half2_"+gene+"_seqs"] = flatten(half2_gene_seq_holder)
    globals()[jumpfunction] = flatten(jumpfunction_holder)

    # Build Aho-Corasick tries
    globals()[gene+"_builder"] = AcoraBuilder()
    for i in range(0,len(globals()[gene+"_seqs"])):
        globals()[gene+"_builder"].add(str(globals()[gene+"_seqs"][i])) # Add all V tags to keyword trie
    globals()[gene+"_key"] = globals()[gene+"_builder"].build()

    # And tries for split, half-tags
    globals()[gene+"_half1_builder"] = AcoraBuilder()
    for i in range(0,len(globals()["half1_"+gene+"_seqs"])):
        globals()[gene+"_half1_builder"].add(str(globals()["half1_"+gene+"_seqs"][i]))
    globals()["half1_"+gene+"_key"] = globals()[gene+"_half1_builder"].build()

    globals()[gene+"_half2_builder"] = AcoraBuilder()
    for i in range(0,len(globals()["half2_"+gene+"_seqs"])):
        globals()[gene+"_half2_builder"].add(str(globals()["half2_"+gene+"_seqs"][i]))
    globals()["half2_"+gene+"_key"] = globals()[gene+"_half2_builder"].build()

  return chain_order

def get_v_deletions( read, v_match, temp_end_v, v_regions_cut ):
    # This function determines the number of V deletions in sequence read
    # by comparing it to v_match, beginning by making comparisons at the
    # end of v_match and at position temp_end_v in read.
    function_temp_end_v = temp_end_v
    pos = len(v_regions_cut[v_match]) -10    # changed from -1 for new checking technique
    is_v_match = 0
    
    # Catch situations in which the temporary end of the V exists beyond the end of the read
    if function_temp_end_v >= len(read):
      counts['v_del_failed_tag_at_end'] += 1
      return
    
    function_temp_end_v += 1
    num_del = 0

    while is_v_match == 0 and 0 <= function_temp_end_v < len(read):
        # Require a 10 base match to determine where end of germ-line sequence lies
        if str(v_regions_cut[v_match])[pos:pos+10] == read[function_temp_end_v-10:function_temp_end_v]:
            is_v_match = 1
            deletions_v = num_del            
            end_v = temp_end_v - num_del
        else:
            pos -= 1
            num_del += 1
            function_temp_end_v -= 1

    if is_v_match == 1:
        return [end_v, deletions_v]
    else:
        counts['v_del_failed'] += 1
        return 

def get_j_deletions( read, j_match, temp_start_j, j_regions_cut ):
    # This function determines the number of J deletions in sequence read
    # by comparing it to j_match, beginning by making comparisons at the
    # end of j_match and at position temp_end_j in read.
    function_temp_start_j = temp_start_j
    pos = 0
    is_j_match = 0
    while is_j_match == 0 and 0 <= function_temp_start_j+2 < len(str(read)):
        # Require a 10 base match to determine where end of germ-line sequence lies
        if str(j_regions_cut[j_match])[pos:pos+10] == read[function_temp_start_j:function_temp_start_j+10]:
            is_j_match = 1
            deletions_j = pos
            start_j = function_temp_start_j
        else:
            pos += 1
            function_temp_start_j += 1
            
    if is_j_match == 1:
        return [start_j, deletions_j]
    else:
        counts['j_del_failed'] += 1
        return 

def get_v_tags(file_v, half_split):
    #"""Read V tags in from file"""
    v_seqs = [] # Holds all V tags
    jump_to_end_v = [] # Holds the number of jumps to make to look for deletions for each V region once the corresponding tag has been found
    for line in file_v:
        elements = line.rstrip("\n") # Turns every element in a text file line separated by a space into elements in a list
        v_seqs.append(string.split(elements)[0]) # Adds elements in first column iteratively
        jump_to_end_v.append(int(string.split(elements)[1])) # Adds elements in second column iteratively

    half1_v_seqs = []
    half2_v_seqs = []

    for i in range(len(v_seqs)):
        half1_v_seqs.append(v_seqs[i][0:half_split])
        half2_v_seqs.append(v_seqs[i][half_split:])
    
    return [v_seqs, half1_v_seqs, half2_v_seqs, jump_to_end_v]

def get_j_tags(file_j, half_split):
    """Read J tags in from file"""
    j_seqs = [] # Holds all J tags
    jump_to_start_j = [] # Holds the number of jumps to make to look for deletions for each J region once the corresponding tag has been found

    for line in file_j:
        elements = line.rstrip("\n")
        j_seqs.append(string.split(elements)[0])
        jump_to_start_j.append(int(string.split(elements)[1]))

    half1_j_seqs = []
    half2_j_seqs = []

    for j in range(len(j_seqs)):
        half1_j_seqs.append(j_seqs[j][0:half_split])
        half2_j_seqs.append(j_seqs[j][half_split:])

    return [j_seqs, half1_j_seqs, half2_j_seqs, jump_to_start_j]

def findTCRs(fqfile, write_type):
    # Scroll through input file and find TCRs
  print "Writing to "+name_results+suffix+"..."
  with open(name_results + suffix, write_type) as outfile:   
    with opener(fqfile) as f:
      
      for readid, seq, qual in readfq(f):
        start_time = time()
        
        if inputargs['nobarcoding'] == False:
          bc = seq[:30]   
          vdj = seq[30:] 
        else:
          vdj = seq

        if inputargs['nobarcoding'] == False:
          if "N" in bc and inputargs['allowNs'] == False:       # Ambiguous base in barcode region
            counts['dcrfilter_barcodeN'] += 1
        
        counts['read_count'] += 1
        if counts['read_count'] % 100000 == 0 and inputargs['dontcount'] == False:
          print '\t read', counts['read_count'] 
    
        # Get details of the VJ recombination

        if inputargs['orientation'] == 'reverse':
          frameR = 'reverse'
          recomR = dcr(revcomp(vdj), inputargs, chain_order)
          recomF = None

        elif inputargs['orientation'] == 'forward':
          frameF = 'forward'
          recomF = dcr(vdj, inputargs, chain_order)
          recomR = None

        elif inputargs['orientation'] == 'either':              # Looks for reverse, but will look for forward if no reverse found
          recomR = dcr(revcomp(vdj), inputargs, chain_order)
          frameR = 'reverse'
          recomF = None
          if not recomR:
            recomF = dcr(vdj, inputargs, chain_order)
            frameF = 'forward'
            recomR = None

        elif inputargs['orientation'] == 'both':
          recomR = dcr(revcomp(vdj), inputargs, chain_order)
          frameR = 'reverse'
          recomF = dcr(vdj, inputargs, chain_order)
          frameF = 'forward'

        if recomR:
          counts['vj_count'] += 1
          dcr_string = build_dcr_string(recomR, frameR, qual, readid, stemplate)
          outfile.write(dcr_string + '\n')
       
        if recomF:        
          counts['vj_count'] += 1
          dcr_string = build_dcr_string(recomF, frameF, qual, readid, stemplate)
          outfile.write(dcr_string + '\n')



def build_dcr_string(recom, frame, qual, readid, stemplate):

 # vdjqual = qual[30:]

  if frame == 'reverse':                      # note: this handles only nbc case
    tcrQ = qual[::-1][recom[3]:recom[4]]
  elif frame == 'forward':
    tcrQ = qual[recom[3]:recom[4]]

  if inputargs['nobarcoding'] == False:
    bcQ = qual[:30]
    dcr_string = stemplate.substitute( chain = str(recom[5]) + ',', v = str(recom[0]) + ',', j = str(recom[1]) + ',', del_v_or_j = str(recom[2]) + ',', \
    seqid = readid + ',', tcr_seq = str(recom[3]) + ',', \
    tcr_qual = tcrQ + ',', barcode = bc + ',', barqual = bcQ )
    return dcr_string

  else:
    dcr_string = stemplate.substitute(chain = str(recom[5]) + ',', v = str(recom[0]) + ',', j = str(recom[1]) + ',', seqid = readid + ',' , tcr_seq = str(recom[2]) + ',', tcr_qual = tcrQ)   
    return dcr_string


def sort_permissions(fl):
  # Need to ensure proper file permissions on output data
    # If users are running pipeline through Docker might otherwise require root access
  if oct(os.stat(fl).st_mode)[4:] != '666':
    os.chmod(fl, 0o666)


##########################################################
############# READ IN COMMAND LINE ARGUMENTS #############
##########################################################

if __name__ == '__main__':
  s_t = time()

  inputargs = vars(args())
  
  print "Running Decombinator version", __version__

  # Determine compression status (and thus opener required)
  if inputargs['fastq'].endswith('.gz'):
    opener = gzip.open
  else:
    opener = open

  # Brief FASTQ sanity check
  if inputargs['dontcheck'] == False:
    if fastq_check(inputargs['fastq']) <> True:
      print "FASTQ sanity check failed reading", inputargs['fastq'], "- please ensure that this file is a properly formatted FASTQ."
      sys.exit()
  
  # Get TCR gene information
  chain_order = import_tcr_info(inputargs)
  
  counts['start_time'] = time()
  
  #########################################################
  ############# SCROLL THROUGH FILE & ANALYSE #############
  #########################################################

  print "Decombining FASTQ data..."

  suffix = "." + inputargs['extension']

  bnam1 = os.path.basename(inputargs['fastq'])
  snam1 = bnam1.split(".")[0]

  if inputargs['fastq2']: #naming hack if two reads are included in input. Works with currently name data files separated by "_"
    bnam2 = os.path.basename(inputargs['fastq2'])
    snam2 = bnam2.split(".")[0]
    samplenam = "_".join(collections.OrderedDict.fromkeys((snam1+"_"+snam2).split("_")).keys())
  else:
    samplenam = snam1

  # If chain had not been autodetected, write it out into output file
  if counts['chain_detected'] == 1:
    name_results = inputargs['prefix'] + samplenam
  else:
    name_results = inputargs['prefix'] + "_".join(map(chainnams.__getitem__, chain)) + "_" + samplenam

  if inputargs['nobarcoding'] == False:
    stemplate = string.Template('$chain $v $j $del_v_or_j $seqid $tcr_seq $tcr_qual $barcode $barqual')
  else:  
    stemplate = string.Template('$chain $v $j $seqid $tcr_seq $tcr_qual')
    found_tcrs = coll.Counter()

  findTCRs(inputargs['fastq'], 'w')

  if inputargs['fastq2']:
    findTCRs(inputargs['fastq2'], 'a')

  counts['end_time'] = time()
  timetaken = counts['end_time']-counts['start_time']

  if inputargs['dontgzip'] == False:
    print "Compressing Decombinator output file,", name_results + suffix, "..."
    
    with open(name_results + suffix) as infile, gzip.open(name_results + suffix + '.gz', 'wb') as outfile:
        outfile.writelines(infile)
    os.unlink(name_results + suffix)

    outfilenam = name_results + suffix + ".gz"
  else:
    outfilenam = name_results + suffix

  sort_permissions(outfilenam)
  
  ##############################################
  ############# WRITE SUMMARY DATA #############
  ##############################################

  print "Analysed", "{:,}".format(counts['read_count']), "reads, finding", "{:,}".format(counts['vj_count']), ", ".join(map(chainnams.__getitem__, chain)), "VJ rearrangements"
  print "Reading from", inputargs['fastq'] + ", writing to", outfilenam
  print "Took", str(round(timetaken,2)), "seconds"

  # Write data to summary file
  if inputargs['suppresssummary'] == False:
    
    # Check for directory and make summary file
    if not os.path.exists('Logs'):
      os.makedirs('Logs')
    date = strftime("%Y_%m_%d")
    
    # Check for existing date-stamped file
    summaryname = "Logs/" + date + "_" + samplenam + "_Decombinator_Summary.csv"
    if not os.path.exists(summaryname): 
      summaryfile = open(summaryname, "w")
    else:
      # If one exists, start an incremental day stamp
      for i in range(2,10000):
        summaryname = "Logs/" + date + "_" + samplenam + "_Decombinator_Summary" + str(i) + ".csv"
        if not os.path.exists(summaryname): 
          summaryfile = open(summaryname, "w")
          break

    # Generate string to write to summary file 
    summstr = "Property,Value\nDirectory," + os.getcwd() + "\nInputFile," + inputargs['fastq'] + "\nOutputFile," + outfilenam \
      + "\nDateFinished," + date + "\nTimeFinished," + strftime("%H:%M:%S") + "\nTimeTaken(Seconds)," + str(round(timetaken,2)) + "\n\nInputArguments:,\n"
    for s in ['species', 'chain','extension', 'tags', 'dontgzip', 'allowNs', 'orientation', 'lenthreshold']:
      summstr = summstr + s + "," + str(inputargs[s]) + "\n"

    counts['pc_decombined'] = counts['vj_count'] / counts['read_count']

    summstr = summstr + "\nNumberReadsInput," + str(counts['read_count']) + "\nNumberReadsDecombined," + str(counts['vj_count']) + "\nPercentReadsDecombined," + str( round(counts['pc_decombined'], 3))

    # Half tag matching details
    summstr = summstr + "\n\nReadsAssignedUsingHalfTags:,\nV1error," + str(counts['verr1']) \
      + "\nV2error," + str(counts['verr2']) \
      + "\nJ1error," + str(counts['jerr1']) \
      + "\nJ2error," + str(counts['jerr2'])
    
    # Number reads filtered out
    summstr = summstr + "\n\nReadsFilteredOut:,\nAmbiguousBaseCall(DCR)," + str(counts['dcrfilter_intertagN']) \
      + "\nAmbiguousBaseCall(Barcode)," + str(counts['dcrfilter_barcodeN']) \
      + "\nOverlongInterTagSeq," + str(counts['dcrfilter_toolong_intertag']) \
      + "\nImpossibleDeletions," + str(counts['dcrfilter_imposs_deletion']) \
      + "\nOverlappingTagBoundaries," + str(counts['dcrfilter_tag_overlap']) \
        
    ##########################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##########################
    summstr = summstr + "\n\nReadsFailedAssignment:,\nMultipleVtagMatches," + str(counts['multiple_v_matches']) \
      + "\nVTagAtEndRead," + str(counts['v_del_failed_tag_at_end']) \
      + "\nVDeletionsUndetermined," + str(counts['v_del_failed']) \
      + "\nFoundV1HalfTagNotV2," + str(counts['foundv1notv2']) \
      + "\nFoundV2HalfTagNotV1," + str(counts['foundv2notv1']) \
      + "\nNoVDetected," + str(counts['no_vtags_found']) \
      + "\nMultipleJTagMatches," + str(counts['multiple_j_matches']) \
      + "\nJDeletionsUndermined," + str(counts['j_del_failed']) \
      + "\nFoundJ1HalfTagNotJ2," + str(counts['foundj1notj2']) \
      + "\nFoundJ2HalfTagNotJ1," + str(counts['foundj2notj1']) \
      + "\nNoJDetected," + str(counts['no_j_assigned']) 
      #+ "\nVJGeneAssignmentFailed," + str(counts['VJ_assignment_failed'])     
        
    print >> summaryfile, summstr 
    summaryfile.close()
    sort_permissions(summaryname)
  print("--- %s seconds ---" % (time() - s_t))
  sys.exit()