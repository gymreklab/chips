# TODO notes on CTCF:
# - peaks in sim data look much broader than in ENCODE data. why?
# - maybe for CTCF we shouldn't do scale-outliers? some peaks tiny in real data are huge in sim data

# To run workflow, edit params in Snakefile if necessary then run:
snakemake

# To make an ENCODE tdf for comparison
igvtools count -z 5 -w 25 -e 0 /storage/mgymrek/chipmunk_round2/encode/GM12878_CTCF_ENCFF124LLI_ENCFF256QBB/GM12878_CTCF_ENCFF124LLI_ENCFF256QBB.flagged.bam /storage/mgymrek/chipmunk_round2/encode/GM12878_CTCF_ENCFF124LLI_ENCFF256QBB/GM12878_CTCF_ENCFF124LLI_ENCFF256QBB.flagged.tdf /storage/resources/dbase/human/hg19/hg19.fa
