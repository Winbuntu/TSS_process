import HTSeq
import numpy
from matplotlib import pyplot

import subprocess



halfwinwidth = 3000
fragmentsize = 200


def get_maped_reads(sort_bam_name):
    '''
    This function use samtools to count how many mapped reads in the bam file
    :param sort_bam_name:
    :return: numeric, number of mapped reads
    '''

    command_line = '''samtools view -F 4''' + sort_bam_name + '''| wc -l'''
    return int(subprocess.check_output(command_line, shell=True))



def construct_gene_Dic(gtf_file_name = "Homo_sapiens.GRCh37.56_chrom1.gtf"):

    '''
    This function construction dict Dic[gene_name] = [p.chrom, TSS_start,TSS_end]
    for construct window.

    :param gtf_file_name:
    :return: Dictionary
    '''
    gtffile = HTSeq.GFF_Reader(gtf_file_name)

    gene_Dic = {}

    for feature in gtffile:
        if feature.type == "exon" and feature.attr["exon_number"] == "1":
            #print str(feature.attr["gene_name"])  # this line get gene name

            if gene_Dic.has_key(feature.attr["gene_name"]) == False:
                gene_Dic[feature.attr["gene_name"]] = set([feature.iv.start_d_as_pos])
            else:
                gene_Dic[feature.attr["gene_name"]].add(feature.iv.start_d_as_pos)


    return gene_Dic


def get_profile_using_gene_Dic(gene_Dic):

    sortedbamfile = HTSeq.BAM_Reader("SRR001432_head_sorted.bam")

    profile = numpy.zeros(2 * halfwinwidth, dtype='i')

    for gene_name in gene_Dic:

        # if gene_name != "WASH5P": continue # select WASH5P gene for fun

        for TSS in gene_Dic[gene_name]:

            window = HTSeq.GenomicInterval(TSS.chrom,
                                           TSS.pos - halfwinwidth - fragmentsize, TSS.pos + halfwinwidth + fragmentsize, "." )

            for almnt in sortedbamfile[window]:
                almnt.iv.length = fragmentsize
                if TSS.strand == "+":
                    start_in_window = almnt.iv.start - TSS.pos + halfwinwidth
                    end_in_window = almnt.iv.end - TSS.pos + halfwinwidth
                else:
                    start_in_window = TSS.pos + halfwinwidth - almnt.iv.end
                    end_in_window = TSS.pos + halfwinwidth - almnt.iv.start
                start_in_window = max(start_in_window, 0)
                end_in_window = min(end_in_window, 2 * halfwinwidth)
                if start_in_window >= 2 * halfwinwidth or end_in_window < 0:
                    continue
                profile[start_in_window: end_in_window] += 1

    pyplot.plot( numpy.arange( -halfwinwidth, halfwinwidth ), profile )
    pyplot.show()







if __name__ == "__main__":
    d = construct_gene_Dic()
    get_profile_using_gene_Dic(d)


#print sys.maxint

#"Homo_sapiens.GRCh37.56_chrom1.gtf"
#
# sortedbamfile = HTSeq.BAM_Reader( "SRR001432_head_sorted.bam" )
# gtffile = HTSeq.GFF_Reader( "Homo_sapiens.GRCh37.56_chrom1.gtf" )
#
#
# tsspos = set()
#
# for feature in gtffile:
#     if feature.type == "exon" and feature.attr["exon_number"] == "1":
#
#         print feature.attr["gene_name"] # this line get gene name
#
#         tsspos.add( feature.iv.start_d_as_pos )
#         print feature.iv.start_d_as_pos
#
# profile = numpy.zeros( 2*halfwinwidth, dtype='i' )
#
# for p in tsspos:
#     window = HTSeq.GenomicInterval( p.chrom,
#          p.pos - halfwinwidth - fragmentsize, p.pos + halfwinwidth + fragmentsize, "." )
#
#
    # for almnt in sortedbamfile[ window ]:
    #     almnt.iv.length = fragmentsize
    #     if p.strand == "+":
    #         start_in_window = almnt.iv.start - p.pos + halfwinwidth
    #         end_in_window	= almnt.iv.end	- p.pos + halfwinwidth
    #     else:
    #         start_in_window = p.pos + halfwinwidth - almnt.iv.end
    #         end_in_window	= p.pos + halfwinwidth - almnt.iv.start
    #     start_in_window = max( start_in_window, 0 )
    #     end_in_window = min( end_in_window, 2*halfwinwidth )
    #     if start_in_window >= 2*halfwinwidth or end_in_window < 0:
    #         continue
    #     profile[ start_in_window : end_in_window ] += 1
#
#
# pyplot.plot( numpy.arange( -halfwinwidth, halfwinwidth ), profile )
# pyplot.show()