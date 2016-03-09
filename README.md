# recurse_app
The program orf_finder.py can be used to read in a multi-sequence fasta file.  It will generate statistics on this file, including the number of sequences, the longest and shortest sequences, and the lengths of these sequences.  It will also search for all valid open reading frames (orfs) in all the sequences, and report the sequence name of the sequence with the longest valid orf.  It will then take the longest valid orf for every sequence, translate the DNA sequence into protein sequence, and output a fasta file with the protein sequence of the longest valid orf for each input sequence (if one is present).   
