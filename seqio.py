from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator


# Sequence codes:
#   R = G/A (purine)
#   Y = T/C (pyrimidine)
#   K = G/T (ketone)
#   M = A/C (amino group)
#   S = G/C (strong interaction)
#   W = A/T (weak interaction)
#   B = G/T/C (not A)
#   D = G/A/T (not C)
#   H = A/C/T (not G)
#   V = G/C/A (not T/U)
# Extended list of sequence characters copied from Sircah:
COMPLEMENT = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N', 'Y':'R',
              'R':'Y', 'U':'A', 'S':'S', 'W':'W', 'K':'M', 'M':'K',
              'B':'V', 'D':'H', 'H':'D', 'V':'B',
              'a':'t', 'c':'g', 't':'a', 'g':'c', 'n':'n', 'y':'r',
              'r':'y', 'u':'a', 's':'s', 'w':'w', 'k':'m', 'm':'k',
              'b':'v', 'd':'h', 'h':'d', 'v':'b' }

def fasta_read(fastafile):
    for record in SeqIO.parse(fastafile, "fasta"):
        yield record

def fastq_read(fastqfile):
    for title,seq,qual in FastqGeneralIterator(fastqfile):
        yield title,seq,qual

def revseq(string):
    tmp = []
    for each in string:
        tmp.append(each)
    tmp.reverse()
    return "".join(tmp)


def rev_complement(string):
    roc_str=''.join([COMPLEMENT[each] for each in string][::-1])
    return roc_str


def complement(string):
    com_str=''.join([COMPLEMENT[each] for each in string])
    return com_str
