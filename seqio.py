from Bio import SeqIO
from Bio.Seq import Seq

#import pysam


def fasta_read(fastafile):
	for record in SeqIO.parse(fastafile, "fasta"):
		print dir(record)
		yield record


def revseq(string):
	tmp = []
	for each in string:
		tmp.append(each)
	tmp.reverse()
	return "".join(tmp)


def rev_complement(string):
	rule={'A':'T','T':'A','C':'G','G':'C','N':'N','n':'N','a':'T','t':'A','c':'G','g':'C'}
	roc_str=''.join([rule[each] for each in string][::-1])
	return roc_str


def complement(string):
	rule={'A':'T','T':'A','C':'G','G':'C','N':'N','n':'N','a':'T','t':'A','c':'G','g':'C'}
	com_str=''.join([rule[each] for each in string])
	return com_str
