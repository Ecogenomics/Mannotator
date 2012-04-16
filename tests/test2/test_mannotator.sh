#! /bin/bash

#module load blast/2.2.22

CONTIGS_FASTA=contigs.fa
PRODIGAL_GFF=prodigal.gff
RAST_GFF=rast.gff
THREADS=1

BINDIR="../../bin/"
PROG="mannotator"

# Mannotator using nr
echo "Running Mannotator with the nr database..."
PROTDB=./db/Nr_truncated
I2A=./db/Nr_mappings_truncated.txt
OUT=mannotatored_nr.gff
CMD="${BINDIR}${PROG} -n -gffs ${RAST_GFF},${PRODIGAL_GFF} -contigs ${CONTIGS_FASTA} -protdb ${PROTDB} -i2a ${I2A} --seq-embed -threads $THREADS -out $OUT -keep"
echo $CMD
eval $CMD
echo ""


# Mannotator using uniref
echo "Running Mannotator with the uniref database..."
PROTDB=./db/Uniref_truncated
I2A=./db/Uniref_mappings_truncated.txt
OUT=mannotatored_uniref.gff
CMD="${BINDIR}${PROG} -n -gffs ${RAST_GFF},${PRODIGAL_GFF} -contigs ${CONTIGS_FASTA} -protdb ${PROTDB} -i2a ${I2A} --seq-embed -threads $THREADS -out $OUT -keep"
echo $CMD
eval $CMD
echo ""

