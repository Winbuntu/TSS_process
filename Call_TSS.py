import HTSeq
import numpy
from matplotlib import pyplot

import subprocess
import sys




def get_maped_reads(sort_bam_name):
    '''
    This function use samtools to count how many mapped reads in the bam file
    :param sort_bam_name:
    :return: numeric, number of mapped reads
    '''

    command_line = '''samtools view -F 4''' + sort_bam_name + '''| wc -l'''
    return int(subprocess.check_output(command_line, shell=True))



def construct_gene_Dic(gtf_file_name):
    '''
    This function construction dict Dic[gene_name] = [p.chrom, TSS_start,TSS_end]
    for construct window.

    :param gtf_file_name:
    :return: Dictionary
    '''
    pass

#print sys.maxint



sortedbamfile = HTSeq.BAM_Reader( "SRR001432_head_sorted.bam" )
gtffile = HTSeq.GFF_Reader( "Homo_sapiens.GRCh37.56_chrom1.gtf" )
halfwinwidth = 3000
fragmentsize = 200

tsspos = set()

for feature in gtffile:
    if feature.type == "exon" and feature.attr["exon_number"] == "1":

        print feature.attr["gene_name"] # this line get gene name

        tsspos.add( feature.iv.start_d_as_pos )
        print feature.iv.start_d_as_pos

profile = numpy.zeros( 2*halfwinwidth, dtype='i' )

for p in tsspos:
    window = HTSeq.GenomicInterval( p.chrom,
         p.pos - halfwinwidth - fragmentsize, p.pos + halfwinwidth + fragmentsize, "." )


    for almnt in sortedbamfile[ window ]:
        almnt.iv.length = fragmentsize
        if p.strand == "+":
            start_in_window = almnt.iv.start - p.pos + halfwinwidth
            end_in_window	= almnt.iv.end	- p.pos + halfwinwidth
        else:
            start_in_window = p.pos + halfwinwidth - almnt.iv.end
            end_in_window	= p.pos + halfwinwidth - almnt.iv.start
        start_in_window = max( start_in_window, 0 )
        end_in_window = min( end_in_window, 2*halfwinwidth )
        if start_in_window >= 2*halfwinwidth or end_in_window < 0:
            continue
        profile[ start_in_window : end_in_window ] += 1


pyplot.plot( numpy.arange( -halfwinwidth, halfwinwidth ), profile )
pyplot.show()