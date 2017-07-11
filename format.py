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

class GTF_Line(object) :
    """Encapsulates a single line/record in a GTF file."""
    def __init__(self, rawstring) :
        ignorepos  = rawstring.find(COMMENT_START)
        s = rawstring[:ignorepos] if ignorepos >= 0 else rawstring.strip()
        self.parts = s.split('\t')
        if len(self.parts) != ATTR_INDEX+1 :
            raise ValueError('Bad GTF record: had %d columns where %d are required' % (len(self.parts), ATTR_INDEX+1))
        self.attrs = self.setAttr(CARED_ATTR)
        if not (self.attrs[GENE_NAME_ATTR] and self.attrs[GENE_ID_ATTR]) :
            raise ValueError('Missing gene name and gene id')
        elif not self.attrs[GENE_NAME_ATTR] :
            self.attrs[GENE_NAME_ATTR] = self.attrs[GENE_ID_ATTR]
        elif not self.attrs[GENE_ID_ATTR] :
            self.attrs[GENE_ID_ATTR] = self.attrs[GENE_NAME_ATTR]
                        
        self.parts[SEQNAME_INDEX] = self.parts[SEQNAME_INDEX].lower()
        
     
    def attributes(self) :      return self.attrs
    def seqname(self) :         return self.parts[SEQNAME_INDEX]
    def source(self) :          return self.parts[SOURCE_INDEX]
    def feature(self) :         return self.parts[FEATURE_INDEX]
    def start(self) :           return int(self.parts[START_INDEX])
    def end(self) :             return int(self.parts[END_INDEX])
    def score(self) :           return self.parts[SCORE_INDEX]
    def strand(self) :          return self.parts[STRAND_INDEX]
    def strand(self) :          return self.parts[STRAND_INDEX]
    def frame_index(self) :     return self.parts[FRAME_INDEX]
    def raw_attributes(self) :  return self.parts[ATTR_INDEX]
    
    def gene_id(self) :         return self.attrs[GENE_ID_ATTR]
    def gene_name(self) :       return self.attrs[GENE_NAME_ATTR]
    def gene_type(self) :       return self.getAttr[GENE_BIO_ATTR]
    def transcript_id(self) :   return self.getAttr(TRANS_ID_ATTR)
    def transcript_name(self) : return self.getAttr(TRANS_NAME_ATTR)

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

