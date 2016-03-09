# this program will import a multi-sequence fasta file and output a
# series of statistics about that file as well as a fasta file
# with the longest orf in a given reading frame for all sequences

# first open the file
try:
    file_name = open('sample_sequences.fasta','r')
except IOError:
    print ("Can't find your file")

# now build a dictionary to put the sequences in
seqs = {}
for line in file_name:
    line = line.rstrip() # remove return characters from end of line
    if line[0] == '>':  
        words = line.split()
        name = words[0][1:] # remove '>' from name
        seqs [name] = ''  
    else: 
        seqs[name] = seqs[name] + line
file_name.close()

# now build lists with entries from the dictionary
List_of_sequence_names = list(seqs.keys())
List_of_sequences = list(seqs.values())

# print the # of entries in the fasta file
print('\nThere are %i entries in this fasta file' % len(List_of_sequence_names)) 

# now build list with lengths of all sequences
List_of_lengths = []
for entry in List_of_sequences:
    List_of_lengths.append(len(entry))

#compute max, min, and # with max and min values
max_length = max(List_of_lengths)
no_of_maxes = List_of_lengths.count(max(List_of_lengths))
min_length = min(List_of_lengths)
no_of_mins = List_of_lengths.count(min(List_of_lengths))

# now identify the IDs of sequences with max and min lengths
max_names = []
min_names = []
for seq_name in List_of_sequence_names:
    if len(seqs[seq_name]) == max_length:
        max_names.append(seq_name)
    elif len(seqs[seq_name]) == min_length:
            min_names.append(seq_name)

# print stats on sequence lengths
print ('\nMax sequence length: %s bp, Sequences of this length: %s' % (max_length,", ".join(max_names)))
print ('\nMin sequence length: %s bp, Sequences of this length: %s' % (min_length,", ".join(min_names)))

# now use function GetORFs to identify all valid orfs in our dictionary in with reading frame = 1
def GetORFs(dict_sequences, reading_frame = 1):
    """GetORFs takes a dictionary of sequence names with sequences and returns a dictionary of sequence
    names with a list of all valid orfs found within a given frame search (default = 1)."""
    acceptable_stops = ['taa', 'tag', 'tga'] 
    all_orfs_found = {} # set empty dictionary to put found ORF sequences associated with each name
    for sequence_name in dict_sequences:
        all_orfs_found[sequence_name] = [] # add sequence name to output dictionary with orf list blank
        start_locations = [] # set empty list of start codon locations
        stop_locations = [] # set empty list of stop codon locations   
        running_sequence = [] # set empty list to keep track of codons as i scan through
        sequence_to_scan = dict_sequences[sequence_name] # retrieve the string with the sequence for this ID
        range_to_search = range(reading_frame - 1, len(sequence_to_scan),3) # define where we will search for codons
        within_orf = 0 # initialize orf finder
        for i in range_to_search:
            codon_to_search = sequence_to_scan.lower()[i:i+3] # this is the codon we will check
            if codon_to_search == 'atg': # search for start codon
                within_orf = 1 # indicate that a start codon was found
            if within_orf == 1: # if we are inside an orf, keep track of the sequence we just scanned
                running_sequence.append(codon_to_search) # add current codon to running sequence tally
            if codon_to_search in acceptable_stops and within_orf == 1: # if we find a stop and have already found a start codon
                all_orfs_found[sequence_name].append(''.join(running_sequence))  # add the running sequence to the list of found orfs in the orf dictionary
                running_sequence = [] # clear running sequence since we completed an orf
                within_orf = 0 # indicate that the orf is over 
    return all_orfs_found   


#run GetORFs on our dictionary to get all possible orfs in frame 1
all_seq_orfs = GetORFs(seqs,1)

# find the longest orf for each individual sequence ID and for all sequence IDs
orf_max_length = 0 # initialize variable to keep track of longest orf
longest_orfs_dict = {} # this dictionary will keep track of longest orf for each ID
for seq_key in all_seq_orfs: #for every sequence ID
    ID_max_length = 0 # initialize variable to keep track of longest orf for each ID
    for orf in all_seq_orfs[seq_key]: # loop through all orfs for that sequence
        orf_length = len(orf)
        if orf_length > ID_max_length: # if this is the longest orf for this ID so far
            longest_orfs_dict[seq_key] = orf
            ID_max_length = orf_length # set this length to the local max
        if orf_length > orf_max_length:
            orf_max_length = orf_length # set this length to the longest overall orf so far
            orf_max_ID = seq_key

# translate the dna sequence for the longest orf for each ID to protein sequence
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
longest_proteins= {}
for ID in longest_orfs_dict:
    dna_sequence = Seq(str(longest_orfs_dict[ID]), generic_dna)
    longest_proteins[ID]= str(dna_sequence.translate())

# output the longest orf for each ID into a fasta file called proteins_output
output_file = open("proteins_output.fasta", "w")
for k in longest_proteins:
    output_file.write(">" + k + "\n" + longest_proteins[k] + "\n")
output_file.close()

# print stats on orf lengths
print ('\nThe longest orf is for %s: %s (length = %i amino acids)' % (orf_max_ID,longest_proteins[orf_max_ID],orf_max_length))
print ('\nThe  file proteins_output.fasta was generated with the longest protein for each sequence ID
         
                    
