#! /bin/bash

module load blast/2.2.22

CONTIGS_FASTA=./contigs.fa
RAST_GFF=./rast.gff
PRODIGAL_GFF=./prodigal.gff
PROTDB=./db/Nr_truncated
I2A=./db/Nr_mappings_truncated.txt
THREADS=1

# Only split input files and recombine gff files for each contig
#CMD="mannotator.pl -gffs ${RAST_GFF},${PRODIGAL_GFF} -contigs ${CONTIGS_FASTA} -protdb ${PROTDB} -i2a ${I2A} -seq_embed -threads $THREADS -keep -three -four -five"
CMD="mannotator.pl -gffs ${RAST_GFF},${PRODIGAL_GFF} -contigs ${CONTIGS_FASTA} -protdb ${PROTDB} -i2a ${I2A} -seq_embed -threads $THREADS -keep"

echo $CMD
eval $CMD
