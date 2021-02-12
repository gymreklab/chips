import pysam

#bamfile = "/storage/agoren/nextseq-runs/2017-09-19/aligned/hg19/final-bams/F1796-HNFa3_S21.sorted.flagged.bam"
#bamfile = "/storage/kycho/KC231_PolII_LCN/aligned/F1796-GTX_06_1K_S6.sorted.flagged.bam"
#bamfile = "/storage/kycho/KC231_PolII_LCN/aligned/F1796-GTX_5_1K_Hela9_S8.sorted.flagged.bam"
bamfile = "/storage/pandaman/project/asimon-data/ENCFF431YXJ.marked.bam"
count_filtered = 0
with pysam.AlignmentFile(bamfile, "rb") as samfile:
    ref_names = []
    for rn in samfile.references:
        #if ("_" not in rn) and (rn != "chrM"):
        ref_names.append(rn)

    frag_count_dict = {}
    for rn in ref_names:
        for frag in samfile.fetch(rn, 0, samfile.get_reference_length(rn)):
            if (frag.is_paired) and\
                    (frag.is_proper_pair) and\
                    (not frag.is_unmapped) and\
                    (not frag.mate_is_unmapped):
                #if (not frag.is_unmapped) :
                #    count_filtered += 1
                if (not frag.is_read1) and (not frag.is_read2):
                    print (frag)
                    exit(1)
                frag_ref = frag.reference_name
                frag_start = frag.reference_start
                read_len = frag.reference_length
                frag_mate_start = frag.next_reference_start
                read_len_infer = frag.infer_read_length()
                frag_len = frag.template_length
                print (frag_ref, frag_start, read_len, read_len_infer, frag_len, frag_mate_start, frag.is_read1, frag.is_read2)
               
'''
                #frag_key = (frag_ref,\
                #                min(frag_start, frag_mate_start),\
                #                max(frag_start, frag_mate_start))
                frag_key = (frag_ref, frag_start)
                if frag_key in frag_count_dict:
                    frag_count_dict[frag_key] += 1
                else:
                    frag_count_dict[frag_key] = 1

count_dict = {}
for key in frag_count_dict:
    count = frag_count_dict[key]
    if count in count_dict:
        count_dict[count] += 1
    else:
        count_dict[count] = 1

print (count_dict)
print ("filtered:"+str(count_filtered))
'''
