import pysam
import sys
import os

def Tabix(tabfile,chrpos=0,start=3,end=4):
    if os.path.exists(tabfile+".tbi"):
        tabix = pysam.TabixFile(tabfile)
    else:
        try:
            pysam.tabix_index(tabfile,force=True, seq_col=chrpos, start_col=start, end_col=end,)    
            tabix = pysam.TabixFile(tabfile+".gz")
        except IOError,e:
            sys.stderr.write("Failed build index,please sort the file based on contig and position first\n")
    return tabix


def mappingstat(inbam,):
	bamfile = pysam.AlignmentFile(inbam,'rb')
	test=pysam.AlignmentFile('xx.bam','wb',template=bamfile)
	total_reads = Porper_reads = MAPQ5 = MAPQ20 = unmapped = uniquemapped = multimapped = Discordantmapped = 0
	for segment in bamfile:
		if segment.is_supplementary:
			continue
		if segment.is_unmapped:
			unmapped+=1
			total_reads+=1
			continue
		if segment.is_secondary:
			multimapped+=1
			test.write(segment)
			continue
		if segment.is_proper_pair:
			Porper_reads+=1
			total_reads+=1
		else:
			Discordantmapped+=1
			total_reads+=1
		if segment.mapping_quality >= 20:
			MAPQ20+=1
			MAPQ5+=1
		elif segment.mapping_quality >=5:
			MAPQ5+=1
		else:
			pass
	bamfile.close()
	test.close()
	uniquemapped = total_reads - multimapped 
	# the multimapped may not exactly the READS NUM in original sequence data,since a read could have more than 2 hits.Yet in this way,a read was counted more than one time when we multimapped+=1, hence the uniquemapped result from the equatation is not exactly right either.But I think this is not a big deal.
	return total_reads,Porper_reads,MAPQ5,MAPQ20,unmapped,multimapped,Discordantmapped


def singlebaseCov(bamfile,region):
	chr,start,end=region.split("-")
	start=int(start)
	end=int(end)
	samfile=pysam.AlignmentFile(bamfile,"rb")
	basecount=[0,]*100000
	for pileupcolumn in samfile.pileup(chr,start,end,truncate=True):
		basecount[pileupcolumn.n]+=1
	samfile.close()
	return basecount

		
