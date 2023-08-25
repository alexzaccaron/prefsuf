
import argparse
from Bio import SeqIO


# function that takes two sequences seq1 and seq2, and merges them if they overlap.
# from https://www.geeksforgeeks.org/python-program-to-check-overlapping-prefix-suffix-in-two-lists/
def get_suffix_prefix(seq1, seq2, minlen):
	
	for i in range(100):                 # at most this bases are chopped from 5' and 3' ends of seq1 and seq2, respectively
		seq1sub = seq1[:len(seq1)-i]     # get seq1 and seq2 chopped bases at their 5' and 3' ends
		seq2sub = seq2[i:]


		# check if they overlap, with identical suffix and prefix
		for char in range(len(seq1sub)):    
			if seq2sub.startswith(seq1sub[char:]) and len(seq1sub[char:]) >= minlen:

				# if they overlap, return the merged sequence
				return(str(i) + '.' + str(len(seq1sub[char:])))

	return('0') # if they don't overlap, return character 0


# related to argparse to get parameters and show help page
parser = argparse.ArgumentParser(description="This Python script finds the largest suffix of a sequence that is identical to the prefix of another sequence. The input is two files containing sequences in fasta format. The script will perform all pairwise comparisons between the sufffixes of the sequences in file1 and the prefixes of the sequences in file2. It will print to terminal the merged sequences that have the minimum size of identical suffix-prefix.")
parser.add_argument("-1", "--file1", help="Fasta file with sequences to compare their suffixes.")
parser.add_argument("-2", "--file2", help="Fasta file with sequences to compare their prefixes.")
parser.add_argument("-o", "--out", help="Optional. File to output the size of identical suffix-prefix between all pairwise sequences.", required=False)
parser.add_argument("-l", "--min_length", help="Minimum length of suffix-prefix intersection. Default = 200", default=200, type=int)


# getting the parameters
args = parser.parse_args()

# Check if no arguments were provided, and if so, print the help message
if not any(vars(args).values()):
	parser.print_help()
	exit()


file1 = args.file1
file2 = args.file2
outfile = args.out
sp_min_len = args.min_length



def main():

	# open file, if optional output parameter set
	if outfile is not None:
		filehandle = open(outfile, "w")

	# pairwise comparisons, for each seq1 and each seq2
	for record1 in SeqIO.parse(file1, "fasta"):
		for record2 in SeqIO.parse(file2, "fasta"):
			
			# get the size of the longest suffix of seq1 identical to the prefix of seq2
			values = get_suffix_prefix(record1.seq.upper(), record2.seq.upper(), sp_min_len)
			if(len(values)) > 1:
				cutpositions,sp_size = values.split('.')
				cutpositions = int(cutpositions)
				sp_size = int(sp_size)
				print(">" + record1.id + "|" + record2.id + "\n" + 
					record1.seq.upper()[:(len(record1.seq)-cutpositions)] + record2.seq.upper()[(cutpositions+sp_size):])

	if outfile is not None:
		filehandle.close()

		

if __name__ == "__main__":
    main()


