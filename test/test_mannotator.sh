#! /bin/bash

# Give a FASTA contig file as first arg, a RAST GFF file as second arg
#CONTIGS_FASTA=$1
#RAST_GFF=$2

module load blast/2.2.22

CONTIGS_FASTA=contigs.fa
PRODIGAL_GFF=contigs_prodigal.gff

PROTDB=/srv/whitlam/bio/db/ncbi/nr
I2A=/srv/whitlam/bio/db/mannotator/Nr201108_mappings.txt
THREADS=14
#THREADS=1

# Now Mannotator for the annotation
CMD="mannotator.pl -gffs ${PRODIGAL_GFF} -contigs ${CONTIGS_FASTA} -protdb ${PROTDB} -i2a ${I2A} -seq_embed -threads $THREADS"

echo $CMD
eval $CMD
