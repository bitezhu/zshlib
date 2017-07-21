from utils import ezOpen
from utils import sortArr


################################################################################################################################

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
    def frame_index(self):     return self.parts[FRAME_INDEX]
    #If the feature is a coding exon, frame should be a number between 0-2 that represents the reading frame of the first base. If the feature is not a coding exon, the value should be ".".
    def raw_attributes(self):  return self.parts[ATTR_INDEX]
    
    def gene_id(self):         return self.attrs[GENE_ID_ATTR]
    def gene_name(self):       return self.attrs[GENE_NAME_ATTR]
    def gene_type(self):       return self.getAttr(GENE_BIO_ATTR)
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
# BED format  (0-based)
#track name=pairedReads description="Clone Paired Reads" useScore=1
#chr22 1000 5000 cloneA 960 + 1000 5000 0 2 567,488, 0,3512
#chr22 2000 6000 cloneB 900 - 2000 6000 0 2 433,399, 0,3601
#If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray).
#An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line.
#blockSizes - A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
#blockStarts - A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount.


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


HEADER_PFXT      = 'track'
HEADER_PFXB      = 'browser'
ALL_COLUMNS      = [CHROM,CHR_START,CHR_END,NAME,SCORE,STRAND,THICK_START,THICK_END,ITEM_RGB_VALS,BLOCKCOUNT,BLOCKSIZES,BLOCKSTARTS]
REQUIRED_COLUMNS = [CHROM, CHR_START, CHR_END]

class BED_Line(object):
    """Encapsulates a single line/record in a bed12 file."""
    def __init__(self, rawstring):
        parts       = rawstring.strip().split('\t')
        self.attrs  = {}
        self.header = rawstring.startswith(HEADER_PFXT)
        if not self.header:
            self.header = rawstring.startswith(HEADER_PFXB)
        if self.header : return
        for col in ALL_COLUMNS:
            try:
                self.attrs[col] = parts[col]
            except IndexError, ie:
                if col in REQUIRED_COLUMNS:
                    raise ie
                self.attrs[col] = None


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

    def blockspos(self):
        sizes  = [int(x) for x in self.attrs[BLOCKSIZES].strip(',').split(',')]
        starts = [int(x) for x in self.attrs[BLOCKSTARTS].strip(',').split(',')]
        self.blocks = []
        for i,start in enumerate(starts):
            self.blocks.append([start+self.startpos(),start+sizes[i]+self.startpos()])
        return self.blocks



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



######################################################################################################################################
#GFF format (1-based)
#
######################################################################################################################################
#table genePredExt (0-based)
'''
"A gene prediction with some additional info."
(
string name;            "Name of gene (usually transcript_id from GTF)"
string chrom;           "Chromosome name"
char[1] strand;         "+ or - for strand"
uint txStart;           "Transcription start position"
uint txEnd;             "Transcription end position"
uint cdsStart;          "Coding region start"
uint cdsEnd;            "Coding region end"
uint exonCount;         "Number of exons"
uint[exonCount] exonStarts; "Exon start positions"
uint[exonCount] exonEnds;   "Exon end positions"
int score;              "Score"
string name2;           "Alternate name (e.g. gene_id from GTF)"
string cdsStartStat;    "enum('none','unk','incmpl','cmpl')"
string cdsEndStat;      "enum('none','unk','incmpl','cmpl')"
lstring exonFrames;     "Exon frame offsets {0,1,2}"
)
'''
(entrzID_INDEX,txName_INDEX,chromosome_INDEX,strand_INDEX,txStart_INDEX,txEnd_INDEX,CDSStart_INDEX,CDSEnd_INDEX,txExonCount_INDEX,txExonsStart_INDEX,txExonsEnd_INDEX,Score_INDEX,GeneName_INDEX,CDSStartStat_INDEX,CDSEndStat_INDEX,exonFrame_INDEX) = range(16)

COMMENT_START='#'

class GenePredExt(object):
    """Encapsulates a single line/record in a genepred extended file."""
    def __init__(self, rawstring):
        ignorepos = rawstring.find(COMMENT_START)
        s = rawstring[:ignorepos] if ignorepos >= 0 else rawstring.strip()
        self.parts = s.split('\t')
        if len(self.parts) != exonFrame_INDEX+1:
            raise ValueError('Bad genepred record: had %d columns where %d are required' % (len(self.parts), ATTR_INDEX+1))
        self.attrs={}
        self.entrzID       = self.parts[entrzID_INDEX]
        self.txName        = self.parts[txName_INDEX]
        self.chromosome    = self.parts[chromosome_INDEX].lower()
        self.strand        = self.parts[strand_INDEX]
        self.txStart       = int(self.parts[txStart_INDEX])
        self.txEnd         = int(self.parts[txEnd_INDEX])
        self.CDSStart      = int(self.parts[CDSStart_INDEX])
        self.CDSEnd        = int(self.parts[CDSEnd_INDEX]) 
        self.txExonCount   = int(self.parts[txExonCount_INDEX])
        self.txExonsStart  = map(int,self.parts[txExonsStart_INDEX].strip(',').split(','))
        self.txExonsEnd    = map(int,self.parts[txExonsEnd_INDEX].strip(',').split(','))
        self.Score         = int(self.parts[Score_INDEX])
        self.GeneName      = self.parts[GeneName_INDEX]
        self.CDSStartStat  = self.parts[CDSStartStat_INDEX]
        self.CDSEndStat    = self.parts[CDSEndStat_INDEX]
        self.exonFrame     = self.parts[exonFrame_INDEX].strip(',').split(',')
        try:
            assert len(self.txExonsStart) == len(self.txExonsEnd) == self.txExonCount
        except AssertionError,e:
            sys.stderr.write('inconsistent line, maybe the exon number in this record is wrong. \n')
        self.CDSEnd        = self.CDSEnd-3 if self.strand == "+" else self.CDSEnd  
        #  And CDS should not include STOP codon as it's not translated into an amino acid. I think this maybe a little difference with gtf format. And importantly UTR does't include start_codon and stop_codon,either. So we introduce tmpStart in function getCDSs to correct this 
        self.CDSStart      = self.CDSStart if self.strand == "+" else self.CDSStart+3
        self.exons         = []
        self.introns       = []
        self.cds           = []
        self.utr           = []
        self.fp_utr        = []
        self.tp_utr        = []
        self.start_codon   = []
        self.stop_codon    = []

    def getExons(self):
        for i in range(self.txExonCount):
            self.exons.append([self.txExonsStart[i],self.txExonsEnd[i]])
    
    def getDonors(self):
        if self.strand == '+':
            self.donors   = [[x,x+2] for x in self.txExonsEnd]
        else:
            self.accepors = [[x-2,x] for x in self.txExonsStart]
        
    def getAcceptors(self):
        if self.strand == '+':
            self.acceptors = [[x-2,x] for x in self.txExonsStart]
        else:
            self.acceptors = [[x,x+2] for x in self.txExonsEnd]

    def getCDSs(self):
        if self.CDSStartStat is 'unk' or self.CDSEndStat is 'unk':
        # in this case, CDSStart is equal to CDSEnd , so we will not try to infer cds boundrary     
            return
        for exonStart,exonEnd in self.exons:
            if exonEnd < self.CDSStart:
                self.utr.append([exonStart,exonEnd])
                continue
            elif exonStart < self.CDSStart < exonEnd:
                tmpStart = self.CDSStart if self.strand=='+' else self.CDSStart -3
                self.utr.append([exonStart,tmpStart]) 
                self.cds.append([self.CDSStart,exonEnd])
            elif self.CDSStart < exonEnd < self.CDSEnd:
                self.cds.append([exonStart,exonEnd])
            elif exonStart < self.CDSEnd < exonEnd:
                tmpEnd = self.CDSEnd+3 if self.strand=='+' else self.CDSEnd
                self.utr.append([tmpEnd,exonEnd])
                self.cds.append([exonStart,self.CDSEnd])
            else:
                self.utr.append([exonStart,exonEnd])
                continue
        self.cds = sortArr(self.cds,0)
        return True

    def inferUTR(self):
        if self.strand == "+":
            for s,e in self.utr:
                if s < self.CDSStart:
                    self.fp_utr.append([s,e])
                else:
                    self.tp_utr.append([s,e])
        else:
            for s,e in self.utr:
                if s < self.CDSStart:
                    self.tp_utr.append([s,e])
                else:
                    self.fp_utr.append([s,e])
        return True

    def inferCodon(self):
        if not self.cds : return
        if self.strand == '+':
            self.start_codon  = [self.CDSStart,self.CDSStart+3]
            self.stop_codon   = [self.CDSEnd,self.CDSEnd+3]
        else:
            self.start_codon  = [self.CDSEnd-3,self.CDSEnd]
            self.stop_codon   = [self.CDSStart-3,self.CDSStart]
        if self.CDSStartStat is 'incmpl':
            self.start_codon  = []
        if self.CDSEndStat is 'incmpl':
            self.stop_codon   = []
        return True

    def getIntrons(self):
        for i in range(self.txExonCount-1):
            self.introns.append([self.txExonsEnd[i],self.txExonsStart[i+1]-1])
        self.introns=sortArr(self.introns,0)
        return True

