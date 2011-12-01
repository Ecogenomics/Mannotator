#! /bin/bash

module load blast/2.2.22

CONTIGS_FASTA=contigs.fa
PRODIGAL_GFF=contigs_prodigal.gff
THREADS=4


# Mannotator using nr
echo "Running Mannotator with the nr database..."
PROTDB=./db/Nr_truncated
I2A=./db/Nr_mappings_truncated.txt
OUT=mannotatored_nr.gff
CMD="mannotator.pl -gffs ${PRODIGAL_GFF} -contigs ${CONTIGS_FASTA} -protdb ${PROTDB} -i2a ${I2A} -seq_embed -threads $THREADS -blastx -out $OUT"
echo $CMD
eval $CMD


# Mannotator using uniref
echo "Running Mannotator with the uniref database..."
PROTDB=./db/Uniref_truncated
I2A=./db/Uniref_mappings_truncated.txt
OUT=mannotatored_uniref.gff
CMD="mannotator.pl -gffs ${PRODIGAL_GFF} -contigs ${CONTIGS_FASTA} -protdb ${PROTDB} -i2a ${I2A} -seq_embed -threads $THREADS -blastx -out $OUT"
echo $CMD
eval $CMD
