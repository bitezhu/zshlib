CDS_TYPE          = 'CDS'
EXON_TYPE         = 'exon'
UTR_TYPE          i= 'UTR'
GENE_TYPE         = 'gene'
MRNA_TYPE         = 'mrna'
TRANSCIRPT_TYPE   = 'transcript'
START_CODON_TYPE  = 'start_codon'
STOP_CODON_TYPE   = 'stop_codon'

CDS_TYPES       = [FP_UTR_TYPE, TP_UTR_TYPE, CDS_TYPE]

class BaseFeature(object):
    def __init__(self, chromosome, start, end, strand, featureType, attr={}) :
        self.chromosome  = chromosome
        start, end       = map(int,[start,end])
        self.start       = min(start,end)
        self.end         = max(start,end)
        self.strand      = strand
        self.parent      = None
        self.featureType = featureType
        self.attrs  = attr

    def acceptor(self):
    """
    Returns the location where the acceptor dimer begins: 2nt
    upstream of the exon start position.  This is 2 before the
    start on the + strand and the exact start on the - strand.
    + Example: AG|CGTATTC
    - Example: GAATACG|CT (reverse-complement)
    """
        return [self.start-2,self.start] if self.strand == '+' else [self.end,self.end+2]

    def donor(self):
        return [self.end,self.end+2] if self.strand == '+' else [self.start-2,self.start]

    def contains(self, chromosome, pos, strand):
        return (strand == self.strand and chromosome == self.chromosome and self.start <= pos <= self.end)

    def updateAttr(self, altAttributes={}) :
        attrs = dict(self.attrs)
        attrs.update(altAttributes)
        self.attrs = attrs


class Exon(BaseFeature):
    def __init__(self, chromosome, start, end, strand, transcriptid, attr={}):
        BaseFeature.__init__(self, chromosome, start, end, strand, EXON_TYPE, attr)
        self.parent = ''
    
    def _addParent(self, transcriptid):
        self.parent = transcriptid
    

class CDS(Exon):
    def __init__(self, chromosome, start, end, strand, transcriptid, attr={}):
        BaseFeature.__init__(self, chromosome, start, end, strand, CDS_TYPE, attr)
        self.parent = ''

    def _addParent(self, transcriptid):
        self.parent = transcriptid

class UTR(BaseFeature):
    def __init__(self, chromosome, start, end, strand, transcriptid, attr={}):
        BaseFeature.__init__(self, chromosome, start, end, strand, UTR_TYPE, attr)
        self.parent = ''

    def _addParent(self, transcriptid):
        self.parent = transcriptid

class StartCodon(Exon):
    def __init__(self, chromosome, start, end, strand, transcriptid, attr={}):
        BaseFeature.__init__(self, chromosome, start, end, strand, START_CODON_TYPE, attr)
        self.parent = ''
        self.SC = 1
        self.EC = 0
    def _addParent(self, transcriptid):
        self.parent = transcriptid

class EndCodon(Exon):
    def __init__(self, chromosome, start, end, strand, transcriptid, attr={}):
        BaseFeature.__init__(self, chromosome, start, end, strand, END_CODON_TYPE, attr)
        self.parent = ''
        self.SC = 0
        self.EC = 1   
    def _addParent(self, transcriptid):
        self.parent = transcriptid



class transcript(BaseFeature):
    """
    An mRNA acts like an isoform in that it is associated with a parent gene
    and contains a number of coding sequences (CDS).
    """
    def __init__(self, chromosome, start, end, strand, featureType, id, attr={}):
        BaseFeature.__init__(self, chromosome, start, end, strand, attr)
        self.id          = id
        self.exons       = []
        self.features    = []
        self.cds         = []
        self.cdsMap      = {}
        self.start_codon = None
        self.stop_codon  = None
        self.fp_utr      = []
        self.tp_utr      = []
        self.utrMap      = {}

    def addCDS(self, cds):
        if cds.strand != self.strand :
            raise Exception("ERROR: strand '%s' of CDS from transcript %s does not match gene strand '%s'" % (cds.strand, cds.parent, self.strand))
        if cds.chromosome != self.chromosome :
            raise Exception("ERROR: chromosome '%s' of CDS from transcript %s does not match gene chromosome '%s'" % (cds.chromosome, cds.parent, self.chromosome))
        cdsTuple = (cds.start,cds.end)
        try :
            ignore = self.cdsMap[cdsTuple]
            return False
        except KeyError :
            self.cdsMap[cdsTuple] = cds
        self.cds.append(list(cdsTuple))
        return True
    
    def addUTR(self,utr):
        if utr.strand != self.strand :
            raise Exception("ERROR: strand '%s' of UTR from transcript %s does not match gene strand '%s'" % (cds.strand, cds.parent, self.strand))
        if utr.chromosome != self.chromosome :
            raise Exception("ERROR: chromosome '%s' of UTR from transcript %s does not match gene chromosome '%s'" % (cds.chromosome, cds.parent, self.chromosome))
        utrTuple = (utr.start,utr.end)
        try :
            ignore = self.utrMap[utrTuple]
            return False
        except KeyError :
            self.utrMap[utrTuple] = utr
        if self.strand == "+":
            self.fp_utr.append(list(utrTuple))
        else:
            self.tp_utr.append(list(utrTuple))
        return True

    def addCodon(self,codon):
        if codon.SC:
            self.start_codon = [codon.start,codon.end]
        else:
            self.stop_codon = [codon.start,codon.end]

    def findCodons(self) :
        """
        Infers a transcript's start and end codon positions based on
        the relative positions of UTR and CDS records.
        """
        if self.end_codon and self.start_codon : return
        if not self.cds : return
        self.cds.sort(reverse=(self.strand=='-'))
        prev = self.cds[0]
        for c in self.cds[1:] :
                    530             if not self.start_codon and prev.featureType == FP_UTR_TYPE and c.featureType == CDS_TYPE :
                         531                 self.start_codon = (c.minpos,c.minpos+2) if self.strand == '+' else (c.maxpos-2,c.maxpos)
                          532             elif not self.end_codon and prev.featureType == CDS_TYPE and c.featureType == TP_UTR_TYPE :
                               533                 self.end_codon = (prev.maxpos-2,prev.maxpos) if self.strand == '+' else (prev.minpos,prev.minpos+2)
                                534             prev = c



