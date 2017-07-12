from utils import ezOpen


CDS_TYPE        = 'CDS'
START_TYPE      = 'start_codon'
STOP_TYPE       = 'stop_codon'
UTR_TYPE        = 'UTR'
INTER_TYPE      = 'inter'
INTER_CNS_TYPE  = 'inter_CNS'
INTRON_CNS_TYPE = 'intron_CNS'
EXON_TYPE       = 'exon'
GENE_TYPE       = 'gene'

#VALID_TYPES     = [CDS_TYPE, START_TYPE, STOP_TYPE, FPU_TYPE, TPU_TYPE, INTER_TYPE, INTER_CNS_TYPE, INTRON_CNS_TYPE, EXON_TYPE]
#EXTRA_TYPES     = [GENE_TYPE, TRANS_TYPE, UTR_TYPE, SEC_TYPE]

# Types found in EMBL files that we actually care about:
KNOWN_TYPES     = [CDS_TYPE, EXON_TYPE]

GENE_ID_ATTR    = 'gene_id'
TRANS_ID_ATTR   = 'transcript_id'
EXON_NO_ATTR    = 'exon_number'
GENE_NAME_ATTR  = 'gene_name'
PROTEIN_ID_ATTR = 'protein_id'
TRANS_NAME_ATTR = 'transcript_name'
GENE_SRC_ATTR   = 'gene_source'
GENE_BIO_ATTR   = 'gene_biotype'
TRANS_SRC_ATTR  = 'transcript_source'
EXON_ID_ATTR    = 'exon_id'

SYNONYM         = {GENE_NAME_ATTR:GENE_ID_ATTR, TRANS_NAME_ATTR:TRANS_ID_ATTR}
COMMENT_START   = '#'

CARED_ATTR      = [GENE_ID_ATTR,TRANS_ID_ATTR,GENE_NAME_ATTR,TRANS_NAME_ATTR,GENE_BIO_ATTR] 

(SEQNAME_INDEX, SOURCE_INDEX, FEATURE_INDEX, START_INDEX, END_INDEX, SCORE_INDEX, STRAND_INDEX, FRAME_INDEX, ATTR_INDEX, COMMENT_INDEX) = range(10)

class GTF_Line(object):
    """Encapsulates a single line/record in a GTF file."""
    def __init__(self, rawstring):
        ignorepos = rawstring.find(COMMENT_START)
        s = rawstring[:ignorepos] if ignorepos >= 0 else rawstring.strip()
        self.parts = s.split('\t')
        if len(self.parts) != ATTR_INDEX+1:
            raise ValueError('Bad GTF record: had %d columns where %d are required' % (len(self.parts), ATTR_INDEX+1))
        self.attrs = self.setAttr(CARED_ATTR)
        if not (self.attrs[GENE_NAME_ATTR] and self.attrs[GENE_ID_ATTR]):
            raise ValueError('Missing gene name and gene id')
        elif not self.attrs[GENE_NAME_ATTR]:
            self.attrs[GENE_NAME_ATTR] = self.attrs[GENE_ID_ATTR]
        elif not self.attrs[GENE_ID_ATTR]:
            self.attrs[GENE_ID_ATTR] = self.attrs[GENE_NAME_ATTR]
                        
        self.parts[SEQNAME_INDEX] = self.parts[SEQNAME_INDEX].lower()
        
     
    def attributes(self):      return self.attrs
    def seqname(self):         return self.parts[SEQNAME_INDEX]
    def source(self):          return self.parts[SOURCE_INDEX]
    def feature(self):         return self.parts[FEATURE_INDEX]
    def start(self):           return int(self.parts[START_INDEX]) - 1 # (convert 1-based gtf to 0-based coordiante system)
    def end(self):             return int(self.parts[END_INDEX])
    def score(self):           return self.parts[SCORE_INDEX]
    def strand(self):          return self.parts[STRAND_INDEX]
    def strand(self):          return self.parts[STRAND_INDEX]
    def frame_index(self):     return self.parts[FRAME_INDEX]
    def raw_attributes(self):  return self.parts[ATTR_INDEX]
    
    def gene_id(self):         return self.attrs[GENE_ID_ATTR]
    def gene_name(self):       return self.attrs[GENE_NAME_ATTR]
    def gene_type(self):       return self.getAttr[GENE_BIO_ATTR]
    def transcript_id(self):   return self.getAttr(TRANS_ID_ATTR)
    def transcript_name(self): return self.getAttr(TRANS_NAME_ATTR)

    def getAttr(self, key):
        try:
            return self.attrs[key]
        except KeyError:
            pass
                            
    def setAttr(self,attrList):
        self.attrstr =  self.parts[ATTR_INDEX]
        self.attrs = {}
        for attr in attrList:
            try:
                attrval = self.attrstr.split('%s "'%attr)[1].split(";")[0].strip('"')
                self.attrs[attr] = attrval
            except IndexError:
                self.attrs[attr] = None
        return self.attrs

    def __str__(self):
        attrString = '; '.join(['%s "%s"'%(k,self.attrs[k]) for k in sorted(self.attrs.keys())])
        return '%s\t%s\t%s\t%d\t%d\t.\t%s\t.\t%s' % (self.seqname(), self.source(), self.feature(), self.start(), self.end(), self.strand(), attrString)




#####################################################################################################################################################

CHROM            = 0
CHR_START        = 1
CHR_END          = 2
NAME             = 3
SCORE            = 4
STRAND           = 5
THICK_START      = 6
THICK_END        = 7
ITEM_RGB_VALS    = 8
BLOCKCOUNT       = 9
BLOCKSIZES       = 10
BLOCKSTARTS      = 11


HEADER_PFX       = 'track'
ALL_COLUMNS      = [CHROM,CHR_START,CHR_END,NAME,SCORE,STRAND,THICK_START,THICK_END,ITEM_RGB_VALS,BLOCKCOUNT,BLOCKSIZES,BLOCKSTARTS]


class BED_Line(object):
    """Encapsulates a single line/record in a bed12 file."""
    def __init__(self, rawstring):
        parts       = s.strip().split('\t')
        self.attrs  = {}
        self.header = s.startswith(HEADER_PFX)
        if self.header : return
        for col in ALL_COLUMNS:
            try:
                self.attrs[col] = parts[col]
            except IndexError, ie:
                self.attrs[col] = None
        sizes = [int(x) for x in self.attrs[SIZES].split(',')]
        if self.strand() == '+':
            self.us,self.ds = sizes
        else:
            self.ds,self.us = sizes


    def chromosome(self):
        """Returns the chromosome given by the record."""
        return self.attrs[CHROM].lower()
    
    def startpos(self):
        """Convenience method that returns the first position for the record."""
        return int(self.attrs[CHR_START])

    def endpos(self):
        """Convenience method that returns the last position for the record."""
        return int(self.attrs[CHR_END])
    
    def featurename(self):
        """Returns the name given by the record."""
        return self.attrs[NAME]

    def strand(self):
        """Returns the strand given by the record."""
        return self.attrs[STRAND]

    def acceptor(self):
        """Returns the acceptor site inferred by the record."""
        return self.endpos()-self.ds-1 if self.strand() == '+' else self.startpos()+self.ds

    def donor(self):
        """Returns the donor site inferred by the record."""
        return self.startpos()+self.us if self.strand() == '+' else self.endpos()- self.us-1




###########################################################################################################################
class fastaRecord(object):
    def __init__(self,Id,description,seq):
        self.id = Id
        self.description = description
        self.seq = seq

    def __str__(self):
        return '>' + self.id + ' '+ self.description + '\n' + self.sequence + '\n'

def _fastaHeader(header):
    try:
        Id ,description = header[1:].split(" ",1)
    except ValueError,e:
        Id = header[1:]
        description = None
    return Id,description

def _fasta_itr_from_file(filehandle):
    line = filehandle.readline().strip()
    if line[0] != ">":
        raise ValueError('invalid fasta file,first line must start with ">"')
    seq = []
    Id,description = _fastaHeader(line[1:])
    linemark = 1
    for line in filehandle:
        linemark += 1
        line = line.rstrip()
        if len(line) == 0:
            raise ValueError('there is a blank line %d in fasta file'%linemark)
        if line[0] == ">":
            yield fastaRecord(Id,description,''.join(seq))
            seq = []
            Id ,description = _fastaHeader(line[1:])
            continue
        seq.append(line)
    yield fastaRecord(Id,description,''.join(seq))


def _fasta_itr_from_name(filename):
    "Provide an iteration through the fasta records in the file named fname. "
    f = ezOpen(filename)
    for rec in _fasta_itr_from_file(f):
        yield rec

def _fasta_itr(src):
    """
    Provide an iteration through the fasta records in file `src'.
    Here `src' can be either a file object or the name of a file.
    """
    if type(src) == str:
        return _fasta_itr_from_name(src)
    elif type(src) == file :
        return _fasta_itr_from_file(src)
    else:
        raise TypeError


class fasta_itr(object):
    "An iterator through a sequence of fasta records."
    def __init__(self,src):
        self.__itr = _fasta_itr(src)
       
    def next(self):
        return self.__itr.next()
    
    def __iter__(self):
        return self

