from  utils import sortArr


CDS_TYPE          = 'CDS'
EXON_TYPE         = 'exon'
UTR_TYPE          = 'UTR'
GENE_TYPE         = 'gene'
MRNA_TYPE         = 'mrna'
TRANSCIRPT_TYPE   = 'transcript'
START_CODON_TYPE  = 'start_codon'
STOP_CODON_TYPE   = 'stop_codon'

CDS_TYPES       = [UTR_TYPE, CDS_TYPE]

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
        '''
        + strand example: CGTATTC|GT
        - strand example: AC|GAATACG
        '''
        return [self.end,self.end+2] if self.strand == '+' else [self.start-2,self.start]

    def contains(self, chromosome, pos, strand):
        return (strand == self.strand and chromosome == self.chromosome and self.start <= pos <= self.end)

    def updateAttr(self, altAttributes={}):
        attrs = dict(self.attrs)
        attrs.update(altAttributes)
        self.attrs = attrs
    
    def setAttr(self,Attributes={}):
        attrs = dict(self.attrs)
        for item in Attributes:
            attrs[item] = Attributes[item]
        self.attrs = attrs

    def setParent(self, id):
        self.parent = id

    def __len__(self):
        return self.start-self.end

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

class StopCodon(Exon):
    def __init__(self, chromosome, start, end, strand, transcriptid, attr={}):
        BaseFeature.__init__(self, chromosome, start, end, strand, STOP_CODON_TYPE, attr)
        self.parent = ''
        self.SC = 0
        self.EC = 1   
    def _addParent(self, transcriptid):
        self.parent = transcriptid


def rmExtra(arr):
    unique_arr = []
    [unique_arr.append(obj) for obj in arr if obj not in unique_arr]
    return unique_arr

class transcript(BaseFeature):
    """
    An mRNA acts like an isoform in that it is associated with a parent gene
    and contains a number of coding sequences (CDS).
    """
    def __init__(self, chromosome, start, end, strand, feature, id, attr={}):
        BaseFeature.__init__(self, chromosome, start, end, strand, feature, attr)
        self.id          = id
        self.exons       = []
        self.biotype     = attr['biotype']
        self.feature     = feature
        self.cds         = []
        self.cdsMap      = {}
        self.start_codon = None
        self.stop_codon  = None
        self.utr         = []
        self.fp_utr      = []
        self.tp_utr      = []
        self.utrMap      = {}
        self.exonMap     = {}
        self.introns     = []

    def addexon(self,exon):
        if exon.strand != self.strand:
            raise Exception("ERROR: strand '%s' of exon from transcript %s does not match gene strand '%s'" % (exon.strand, exon.parent, self.strand))
        if exon.chromosome != self.chromosome:
            raise Exception("ERROR: chromosome '%s' of exon from transcript %s does not match gene chromosome '%s'" % (exon.chromosome, exon.parent, self.chromosome))
        exonTuple = (exon.start,exon.end)
        try :
            ignore = self.exonMap[exonTuple]
            return False
        except KeyError,e:
            self.exonMap[exonTuple] = exon
        self.exons.append(list(exonTuple))
        return True

    def addCDS(self, cds):
        if cds.strand != self.strand:
            raise Exception("ERROR: strand '%s' of CDS from transcript %s does not match gene strand '%s'" % (cds.strand, cds.parent, self.strand))
        if cds.chromosome != self.chromosome:
            raise Exception("ERROR: chromosome '%s' of CDS from transcript %s does not match gene chromosome '%s'" % (cds.chromosome, cds.parent, self.chromosome))
        cdsTuple = (cds.start,cds.end)
        try :
            ignore = self.cdsMap[cdsTuple]
            return False
        except KeyError,e:
            self.cdsMap[cdsTuple] = cds
        self.cds.append(list(cdsTuple))
        return True
    
    def addUTR(self,utr):
        if utr.strand != self.strand:
            raise Exception("ERROR: strand '%s' of UTR from transcript %s does not match gene strand '%s'" % (cds.strand, cds.parent, self.strand))
        if utr.chromosome != self.chromosome :
            raise Exception("ERROR: chromosome '%s' of UTR from transcript %s does not match gene chromosome '%s'" % (cds.chromosome, cds.parent, self.chromosome))
        utrTuple = (utr.start,utr.end)
        try :
            ignore = self.utrMap[utrTuple]
            return False
        except KeyError,e:
            self.utrMap[utrTuple] = utr
        self.utr.append(list(utrTuple))
        return True

    def addCodon(self,codon):
        if codon.SC:
            self.start_codon = [codon.start,codon.end]
        else:
            self.stop_codon = [codon.start,codon.end]
        return True

    def inferUTRort(self):
        '''
        Infers a UTR record based on strand and coordinates compare between UTR and cds
        '''
        if not self.utr: return
        self.cds = sortArr(self.cds,0)
        firstcdsStart, firstcdsEnd = self.cds[0]
        for utrs,utre in self.utr:
            if self.strand == "+":
                if utrs < firstcdsStart:
                    self.fp_utr.append([utrs,utre])
                else:
                    self.tp_utr.append([utrs,utre])
            else:
                if utrs < firstcdsStart:
                    self.tp_utr.append([utrs,utre])
                else:
                    self.fp_utr.append([utrs,utre])
        return True

    
    def inferUTR(self):
        '''
        This method would try to infer FP_UTR and TP_UTR
        even without UTR annotation for a transcript. Exons and CDSs are suffcient.
        '''
        if self.fp_utr and self.tp_utr: return
        if not self.cds: return
        self.cds = sortArr(self.cds,0)
        self.exons = sortArr(self.exons,0)
        CDSstart,CDSend = (self.cds[0][0],self.cds[-1][-1]) if self.strand == "+" else (self.cds[-1][-1],self.cds[0][0]) 
        for tmps,tmpe in self.exons:
            if self.strand == "+":
                if tmpe < CDSstart:
                    self.fp_utr.append([tmps,tmpe])
                elif tmps < CDSstart < tmpe:
                    self.fp_utr.append([tmps,CDSstart])
                elif tmps > CDSend:
                    self.tp_utr.append([tmps,tmpe])
                elif tmps < CDSend < tmpe:
                    self.tp_utr.append([CDSend,tmpe])
                else:
                    raise Exception("ERROR: could not determine this exon [%d,%d] is a UTR or not" % (tmps,tmpe,))
            else:
                if tmps > CDSstart:
                    self.fp_utr.append([tmps,tmpe])
                elif tmps < CDSstart < tmpe:
                    self.fp_utr.append([CDSstart,tmpe])
                elif tmpe < CDSend:
                    self.tp_utr.append([tmps,tmpe])
                elif tmps < CDSend < tmpe:
                    self.tp_utr.append([tmps,CDSend])
        self.fp_utr = rmExtra(self.fp_utr)
        self.tp_utr = rmExtra(self.tp_utr)
        return True
    

    def inferCodons(self):
        '''
        This method would try to infer start codon and end codon 
        even without UTR annotation for a transcript.Exons and CDSs are suffcient.
        '''
        if self.stop_codon and self.start_codon: return
        if not self.cds: return
        if self.strand == "+":
            self.start_codon = [self.cds[0][0],self.cds[0][0]+3]
            self.stop_codon = [self.cds[-1][-1],self.cds[-1][-1]+3]
        else:
            self.start_codon = [self.cds[-1][-1]-3,self.cds[-1][-1]]
            self.stop_codon = [self.cds[0][0]-3,self.cds[0][0]]
        return True
    
    def translateRegion(self):
        Start = End = 0
        if self.cds:
            Start = self.cds[0][0]
            End   = self.cds[-1][-1]
        elif self.exons:
            Start = self.exons[0][0]
            End   = self.exons[-1][-1]
        else:
            return
        self.translateRegion = [Start,End]
        return True
    
    def inferIntron(self):
        '''Returns a list of duples containing start/end positions of introns in this transcript.'''
        self.introns = []
        if self.exons:
            intervals = self.exons
        elif self.cds:
            intervals = self.cds
        else:
            return
        for i in xrange(1,len(intervals)) :
            intron = [exons[i-1].end(), exons[i].start()]
            self.introns.append(intron)
        self.introns = sortArr(self.introns,0) 
        return True
    
    def inferRegionLoc(self,region):
        '''Infer a chr:start:end format region's overlap status, Flag '1' mean exonic ,'-1' means intronic, '0' means maybe it's a intergenic region '''
        querychr,querystart,queryend = region.split(":")
        querystart,queryend = map(int,[querystart,queryend])
        exonFlag=0
        for i,(s,e) in enumerate(self.exons):
            if querystart < s < queryend or (querystart>s and queryend<e) or querystart < e < queryend:
                exonFlag = 1
                break
        if not exonFlag:
            for i,(s,e) in enumerate(self.introns):
                if querystart < s < queryend or (querystart>s and queryend<e) or querystart < e < queryend:
                    exonFlag = -1
                    break
        return exonFlag
    
    def upstream(self,distance=1000):
        '''Returns the upstream and downstream flanking region of a transcript record.'''
        uppos   = min(self.cds[0][0],self.exons[0][0])
        downpos = max(self.cds[-1][-1],self.exons[-1][-1])
        if self.strand == '+':
            return [uppos+distance,uppos]
        else:
            return [downpos,downpos+distance]

    def downstream(self,distance=1000):
        uppos   = min(self.cds[0][0],self.exons[0][0])
        downpos = max(self.cds[-1][-1],self.exons[-1][-1])
        if self.strand == '+':
            return [downpos,downpos+distance]
        else:
            return [uppos+distance,uppos]


    def __str__(self):
        return "%s :%s:%d-%d (len=%d, strand=%s), %d exons/cds, translate from %d to %d" % (self.id, self.chromosome, self.start, self.end, self.end-self.start, self.strand, len(self.cds), self.translateRegion[0],self.translateRegion[1])


class Gene(BaseFeature):
    def __init__(self, chromosome, start, end, strand, id, name=None,attr={}):
        BaseFeature.__init__(self, chromosome, start, end, strand, GENE_TYPE , attr)
        self.id          = id
        self.name        = name if name is not None else id
        self.exons       = []
        self.biotype     = attr['biotype']
        self.trancripts  = {}
        self.exons       = []
        self.cds         = []
        self.exonMap     = {}
        self.cdsMap      = {}
        self.string      = ''

        self.start_codons = {}
        self.end_codons   = {}
        
        # All other features associated with genes in an annotation file, such as:
        #    3'/5' UTRs, mRNA, miRNA, siRNA, tRNA, rRNA, ncRNA, snRNA, snoRNA
        self.features = []
    
    def addisoform(self,Newtranscript):
        Newtranscript.setParent(self.id)
        return self.trancripts.setdefault(Newtranscript.id, Newtranscript)
