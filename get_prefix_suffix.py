
import argparse
from Bio import SeqIO


# function that takes two sequences seq1 and seq2, and returns the size of the longest suffix of seq1 that is identical to the prefix of seq2
# from https://www.geeksforgeeks.org/python-program-to-check-overlapping-prefix-suffix-in-two-lists/
def suffix_prefix_size(seq1, seq2):
	for char in range(len(seq1)):           # for each character in seq1
		
		if seq2.startswith(seq1[char:]):    # check if seq2 starts with the suffix of seq1
			res = seq1[char:]               # if yes, return the size of the identical suffix/prefix
			return(len(res))

	return(0)                               # if none found, return 0


# related to argparse to get parameters and show help page
parser = argparse.ArgumentParser(description="This Python script finds the largest suffix of a sequence that is identical to the prefix of another sequence. The input is two files containing sequences in fasta format. The script will perform all pairwise comparisons between the sufffixes of the sequences in file1 and the prefixes of the sequences in file2. It will print to terminal the merged sequences that have the minimum size of identical suffix-prefix.")
parser.add_argument("-1", "--file1", help="Fasta file with sequences to compare their suffixes.")
parser.add_argument("-2", "--file2", help="Fasta file with sequences to compare their prefixes.")
parser.add_argument("-o", "--out", help="Optional. File to output the size of identical suffix-prefix between all pairwise sequences.", required=False)
parser.add_argument("-l", "--min_length", help="Minimum length of suffix-prefix intersection. Default = 200", default=200, type=int)


# getting the parameters
args = parser.parse_args()
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
			sp_size = suffix_prefix_size(record1.seq.upper(), record2.seq.upper())
			
			# print to terminal if size of identical suffix/prefix is at least a minimum
			if sp_size >= sp_min_len:
				print(">" + record1.id + "|" + record2.id + "\n" + 
					record1.seq.upper()[:len(record1.seq)-sp_size] + record2.seq.upper())

			# write size of identical suffix/prefix to file, if optional output parameter set
			if outfile is not None:
				filehandle.write(record1.id + "\t" + record2.id + "\t" + str(sp_size) + "\n")

	# close file, if optional output parameter set
	if outfile is not None:
		filehandle.close()

		

if __name__ == "__main__":
    main()


